# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from .py23_compat import range, int

import numpy as np
from astropy import units as u

from astropy import wcs
from astropy.io import fits
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import Table
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from astropy_healpix import HEALPix
from astropy_healpix import lonlat_to_healpix
from astropy_healpix.healpy import nside2npix

from .interval_set import IntervalSet
from .utils import uniq2orderipix, trailing_zeros

from ...core import PixCoord, PixelRegion, SkyRegion, BoundingBox, RegionMask
from ...core.attributes import RegionMeta, RegionVisual
from ..._geometry.pnpoly import points_in_polygons
from ..._geometry import polygonal_overlap_grid

__all__ = [
    'MOCSkyRegion',
    'MOCPixelRegion',
]

class MOCPixelRegion(PixelRegion):
    """
    A MOC (Multi-Order Coverage map) in pixel coordinates.

    MOC stands for Multi-Order Coverage. It is a spatial and hierarchical description of a region on a sphere.
    MOC is an IVOA standard which was first introduced in the following
    `paper <http://www.ivoa.net/documents/MOCSkyRegion.>`__ .
    MOC are based on the HEALPix sky tessellation using the NESTED numbering scheme. A MOC is a set of
    HEALPix cells at different orders with a maximum resolution corresponding to the order 29 i.e. a cell
    resolution of ~393.2uas.

    * MOCs are usually stored as FITS file containing a list of UNIQ numbers describing the HEALPix cells at \
      different orders. This class aims at creating MOC maps from FITS/json formatted files, FITS images associated with a mask numpy array, \
      an `astropy.coordinates.SkyCoord` object and lon, lat expressed as `astropy.units.Quantity` objects.

    * Basic operations on MOCs are available such as the intersection, union, difference and complement of a MOC.

    * A :meth:`~regions.MOCSkyRegion.contains` method filters positions expressed as \
      `~astropy.units.Quantity` to only keep those lying inside/outside the MOC.

    * A MOC can be serialized to FITS (i.e. a list of UNIQ numbers stored in a \
      binary HDU table) and JSON formats. An optional parameter allows the user to write it to a file.

    Parameters
    ----------
    vertices : `~regions.core.PixCoord` object
        A list of pixel coordinates describing the vertice boundaries of all the HEALPix cells (i.e. 4 pixels for each HEALPix cells).
        Contains two Nx4 numpy arrays, one storing the x coordinates of the HEALPix cell vertices (4 for each cell) and the other for the y coordinates.
    sky_region: `~regions.MOCSkyRegion` object
        The sky region linked to this pixel region object. A MOC is always defined by a MOCSkyRegion first before applying a WCS to it.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    def __init__(self, vertices, sky_region, meta=None, visual=None):
        self.vertices = vertices

        # The backfacing HEALPix cells are culled before plotting the MOC.
        self.vertices_culled = self._backface_culling(self.vertices)
        self.sky_region = sky_region
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('_interval_set',)

    def _backface_culling(self, vertices):
        # Remove cells crossing the MOC after projection
        # The remaining HEALPix cells are used for computing the patch of the MOC
        vx = self.vertices.x
        vy = self.vertices.y

        def cross_product(vx, vy, i):
            cur = i
            prev = (i - 1) % 4
            next = (i + 1) % 4

            # Construct the first vector from A to B
            x1 = vx[:, cur] - vx[:, prev]
            y1 = vy[:, cur] - vy[:, prev]
            z1 = np.zeros(x1.shape)

            v1 = np.vstack((x1, y1, z1)).T
            # Construct the second vector from B to C
            x2 = vx[:, next] - vx[:, cur]
            y2 = vy[:, next] - vy[:, cur]
            z2 = np.zeros(x2.shape)

            v2 = np.vstack((x2, y2, z2)).T
            # Compute the cross product between the two
            return np.cross(v1, v2)

        # A ----- B
        #  \      |
        #   D-----C
        # Compute the cross product between AB and BC
        # and the cross product between BC and CD
        ABC = cross_product(vx, vy, 1)
        CDA = cross_product(vx, vy, 3)

        frontface_cells  = (ABC[:, 2] < 0) & (CDA[:, 2] < 0)

        vx = vx[frontface_cells]
        vy = vy[frontface_cells]

        return vx, vy

    def as_artist(self, origin=(0, 0), **kw_mpl_pathpatch):
        from matplotlib.path import Path
        from matplotlib.patches import PathPatch

        # Only the HEALPix cells facing the camera are used
        vx, vy = self.vertices_culled

        # The border of each HEALPix cells is drawn one at a time
        path_vertices_l = []
        codes = []

        for i in range(vx.shape[0]):
            # Appending to a list is faster that is why we keep
            path_vertices_l += [(vx[i][0], vy[i][0]), (vx[i][1], vy[i][1]), (vx[i][2], vy[i][2]), (vx[i][3], vy[i][3]), (0, 0)]
            codes += [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]

        # Cast to a numpy array
        path_vertices_arr = np.array(path_vertices_l)
        # Add the origin to the path constructed
        path_vertices_arr += np.array(origin)

        path = Path(path_vertices_arr, codes)
        pathpatch = PathPatch(path, **kw_mpl_pathpatch)

        return pathpatch

    def to_sky(self, wcs):
        return self.sky_region

    @property
    def bounding_box(self):
        xmin = self.vertices.x.min()
        xmax = self.vertices.x.max()
        ymin = self.vertices.y.min()
        ymax = self.vertices.y.max()
        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    @property
    def area(self):
        raise NotImplementedError

    def to_mask(self, mode='center', subpixels=5):
        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5
        xmax = float(bbox.ixmax) - 0.5
        ymin = float(bbox.iymin) - 0.5
        ymax = float(bbox.iymax) - 0.5

        vx, vy = self.vertices_culled

        # Loop over all the projeted HEALPix cells to get their overlap grids
        fraction_sum = np.zeros(shape=(ny, nx))
        for i in range(vx.shape[0]):
            fraction = polygonal_overlap_grid(
                xmin, xmax, ymin, ymax,
                nx, ny, vx[i], vy[i],
                use_exact, subpixels,
            )
            # Add the overlap grid of the HEALPix cells to get the
            # overlap area of the total MOC
            fraction_sum += fraction

        # Clip its values to the interval (0, 1)
        fraction_clipped = np.clip(a=fraction_sum, a_min=0, a_max=1)
        return RegionMask(fraction_clipped, bbox=bbox)

    def contains(self, pixcoord):
        """
        Checks whether a position or positions fall inside the region.

        Parameters
        ----------
        pixcoord : `~regions.PixCoord`
            The position or positions to check.
        """
        pixcoord = PixCoord._validate(pixcoord, 'pixcoord')
        x = np.atleast_1d(np.asarray(pixcoord.x, dtype=float))
        y = np.atleast_1d(np.asarray(pixcoord.y, dtype=float))

        shape = x.shape
        vx, vy = self.vertices_culled
        mask = points_in_polygons(x.flatten(), y.flatten(), vx, vy).astype(bool)
        in_poly = mask.reshape(shape)
        if self.meta.get('include', True):
            return in_poly
        else:
            return np.logical_not(in_poly)

class MOCSkyRegion(SkyRegion):
    """
    A MOC (Multi-Order Coverage map)

    MOC stands for Multi-Order Coverage. It is a spatial and hierarchical description of a region on a sphere.
    MOC is an IVOA standard which was first introduced in the following
    `paper <http://www.ivoa.net/documents/MOCSkyRegion.>`__ .
    MOC are based on the HEALPix sky tessellation using the NESTED numbering scheme. A MOC is a set of
    HEALPix cells at different orders with a maximum resolution corresponding to the order 29 i.e. a cell
    resolution of ~393.2uas.

    * MOCs are usually stored as FITS file containing a list of UNIQ numbers describing the HEALPix cells at \
      different orders. This class aims at creating MOC maps from FITS/json formatted files, FITS images associated with a mask numpy array, \
      an `astropy.coordinates.SkyCoord` object and lon, lat expressed as `astropy.units.Quantity` objects.

    * Basic operations on MOCs are available such as the intersection, union, difference and complement of a MOC.

    * A :meth:`~regions.MOCSkyRegion.contains` method filters positions expressed as \
      `~astropy.units.Quantity` to only keep those lying inside/outside the MOC.

    * A MOC can be serialized to FITS (i.e. a list of UNIQ numbers stored in a \
      binary HDU table) and JSON formats. An optional parameter allows the user to write it to a file.

    Parameters
    ----------
    interval_set : `~regions.IntervalSet` object, optional
        A N rows by 2 columns `~numpy.ndarray` storing the set of intervals
        representing the data structure of the MOCSkyRegion.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import MOCSkyRegion

        # Load a MOCSkyRegion from a FITS file including the MOC of the P-GALEXGR6-AIS-FUV survey
        from astropy.utils.data import get_pkg_data_filename
        moc_sky_reg = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits',
                                                                   package='regions'))

        # Configure a WCS
        from regions._utils.examples import make_example_dataset
        config = dict(crpix=(0, 0), crval=(0, 0), cdelt=(-5, 5), shape=(18, 36))
        dataset = make_example_dataset(config=config)

        # Convert the MOCSkyRegion to a MOCPixelRegion for plotting it
        moc_px_reg = moc_sky_reg.to_pixel(dataset.wcs)

        # Define matplotlib to plot the MOCPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        moc_px_reg.plot(ax, edgecolor='red', lw=1)

        plt.xlim(0, 4)
        plt.ylim(4, 8)
        plt.show()
    """
    HPY_MAX_NORDER = 29

    def __init__(self, interval_set=None, meta=None, visual=None):
        interval = IntervalSet() if interval_set is None else interval_set
        self._interval_set = interval
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('_interval_set',)

    def to_pixel(self, wcs):
        healpix_cells_d = self.write(format='json')

        X = np.array([])
        Y = np.array([])
        for order_str, ipix_l in healpix_cells_d.items():
            order = int(order_str)
            # Compute the skycoord boundaries of the HEALPix cells from one order at a time.
            hp = HEALPix(nside=(1 << order), order='nested', frame=ICRS())
            vertices_boundaries = hp.boundaries_skycoord(ipix_l, step=1)

            x, y = skycoord_to_pixel(vertices_boundaries, wcs=wcs)
            if X.size >= 1:
                X = np.vstack((X, x))
                Y = np.vstack((Y, y))
            else:
                X = x
                Y = y

        # Convert them to pixel coordinates
        vertices = PixCoord(X, Y)
        return MOCPixelRegion(vertices=vertices,
                              sky_region=self,
                              meta=self.meta,
                              visual=self.visual)

    def plot(self):
        raise NotImplementedError

    def __eq__(self, another_moc):
        """
        Test equality between self and ``another_moc``

        Parameters
        ----------
        another_moc : `~regions.MOCSkyRegion`
            The moc object to test the equality with.

        Returns
        -------
        result : bool
            True if the interval sets of self and ``another_moc`` are equal.
        """
        if not isinstance(another_moc, MOCSkyRegion):
            raise TypeError('The object you want to test the equality with is not a MOC but a {0}'
                            .format(type(another_moc)))

        return self._interval_set == another_moc._interval_set

    @classmethod
    def from_json(cls, json_moc, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a dictionary of HEALPix cell arrays indexed by their order (i.e. in json format).

        Parameters
        ----------
        json_moc : dict
            A dictionary of (order, [ipix]) key-value pairs representing the set of HEALPix cells defining the MOC map.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object.
        """
        intervals_arr = np.array([], dtype=np.int64)

        for order, pix_l in json_moc.items():
            pix_arr = np.array(pix_l, dtype=np.int64)
            p1 = pix_arr
            p2 = pix_arr + 1
            shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - int(order))

            itv_arr = np.vstack((p1 << shift, p2 << shift)).T
            if intervals_arr.size == 0:
                intervals_arr = itv_arr
            else:
                intervals_arr = np.vstack((intervals_arr, itv_arr))

        return cls(IntervalSet.from_numpy_array(intervals_arr), meta, visual)

    @classmethod
    def from_fits(cls, filename, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a FITS file.

        Works for FITS file in which HEALPix cells are stored as a list of
        UNIQ HEALPix numbers in a binary HDU table. See the `IVOA standard publication part 2.3.1
        <http://www.ivoa.net/documents/MOCSkyRegion.20140602/REC-MOCSkyRegion.1.0-20140602.pdf>`__ for more explanations about NUINQ packing.

        Parameters
        ----------
        filename : str
            The path to FITS file.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object loaded from the FITS file.
        """
        table = Table.read(filename)

        intervals = np.vstack((table['UNIQ'], table['UNIQ']+1)).T

        nuniq_interval_set = IntervalSet.from_numpy_array(intervals)
        interval_set = IntervalSet.from_nuniq_interval_set(nuniq_interval_set)
        return cls(interval_set, meta, visual)

    @classmethod
    def from_image(cls, header, max_norder, mask_arr=None, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from an image stored as a FITS file.

        Parameters
        ----------
        header : `~astropy.io.fits.Header`
            The FITS header of the image.
        max_norder : int
            The maximum order that we want for the new `~regions.MOCSkyRegion`.
        mask_arr : `~numpy.ndarray`, optional
            A 2D boolean numpy array of the same size of the image where pixels evaluated to True are part of
            the new MOC and pixels evaluated to False are not.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object created from the image.
        """
        # load the image data
        height = header['NAXIS2']
        width = header['NAXIS1']

        # use wcs from astropy to locate the image in the world coordinates
        w = wcs.WCS(header)

        if mask_arr is not None:
            # We have an array of pixels that are part of of survey
            y, x = np.where(mask_arr)
            pix_crd = np.dstack((x, y))[0]
        else:
            # If we do not have a mask array we create the moc of all the image
            step_pix = 1
            """
            Coords returned by wcs_pix2world method correspond to pixel centers. We want to retrieve the moc pix
            crossing the borders of the image so we have to add 1/2 to the pixels coords before computing the lonlat.

            The step between two pix_crd is set to `step_pix` but can be diminished to have a better precision at the
            borders so that all the image is covered (a too big step does not retrieve al
            the moc pix crossing the borders of the image).
            """
            x, y = np.mgrid[0.5:(width + 0.5 + step_pix):step_pix, 0.5:(height + 0.5 + step_pix):step_pix]
            pix_crd = np.dstack((x.ravel(), y.ravel()))[0]

        world_pix_crd = w.wcs_pix2world(pix_crd, 1)

        hp = HEALPix(nside=(1 << max_norder), order='nested', frame=ICRS())
        ipix = hp.lonlat_to_healpix(lon=world_pix_crd[:, 0] * u.deg, lat=world_pix_crd[:, 1] * u.deg)
        # remove doubles
        ipix = np.unique(ipix)

        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        # This MOCSkyRegion.will be consistent when one will do operations on the moc (union, inter, ...) or
        # simply write it to a fits or json file
        interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return cls(interval_set, meta, visual)

    @classmethod
    def from_skycoords(cls, skycoords, max_norder, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a `astropy.coordinates.SkyCoord` object.

        Parameters
        ----------
        skycoords : `astropy.coordinates.SkyCoord`
            A set of astropy formatted skycoords.
        max_norder : int
            The maximum order that we want for the new `~regions.MOCSkyRegion`.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object created from a `astropy.coordinates.SkyCoord`
        """
        hp = HEALPix(nside=(1 << max_norder), order='nested')
        ipix = hp.lonlat_to_healpix(skycoords.icrs.ra, skycoords.icrs.dec)

        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return cls(interval_set, meta, visual)

    @classmethod
    def from_lonlat(cls, lon, lat, max_norder, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from ``lon``, ``lat`` `astropy.units.Quantity`.

        Parameters
        ----------
        lon : `astropy.units.Quantity`
            A set of astropy formatted ra quantities.
        lat : `astropy.units.Quantity`
            A set of astropy formatted dec quantities.
        max_norder : int
            The maximum order that we want for the new `~regions.MOCSkyRegion`.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object created from lon, lat `astropy.units.Quantity`.
        """
        hp = HEALPix(nside=(1 << max_norder), order='nested')
        ipix = hp.lonlat_to_healpix(lon, lat)

        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - max_norder)
        intervals_arr = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return cls(interval_set, meta, visual)

    @property
    def max_order(self):
        """
        Returns the deepest order needed to describe self's interval set.
        """
        # TODO: cache value
        combo = int(0)
        for iv in self._interval_set._intervals:
            combo |= iv[0] | iv[1]

        ret = MOCSkyRegion.HPY_MAX_NORDER - (trailing_zeros(combo) // 2)
        if ret < 0:
            ret = 0

        return ret

    @property
    def sky_fraction(self):
        """
        Returns the sky fraction percentage covered by the MOC (between 0 and 1).
        """
        pix_id_arr = self._best_res_pixels()
        nb_pix_filled = pix_id_arr.size
        return nb_pix_filled / float(3 << (2*(self.max_order + 1)))

    def contains(self, ra, dec):
        """
        Check whether (ra, dec) positions fall inside or outside the MOC region.
        The size of ``ra`` must be the same as the size of ``dec``.

        Parameters
        ----------
        ra : `~astropy.units.Quantity`
            A set of astropy formatted ra quantities
        dec: `~astropy.units.Quantity`
            A set of astropy formatted dec quantities

        Returns
        -------
        array : `~numpy.ndarray`
            A boolean numpy array telling which positions are inside/outside of the MOC depending on the
            value of ``self.meta['include']``.
        """
        max_order = self.max_order
        m = np.zeros(nside2npix(1 << max_order), dtype=bool)

        pix_id_arr = self._best_res_pixels()
        m[pix_id_arr] = True

        if self.meta.get('include', False):
            m = np.logical_not(m)

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        pix_arr = hp.lonlat_to_healpix(ra, dec)

        return m[pix_arr]

    def write(self, path=None, format='fits', optional_kw_dict=None, write_to_file=False):
        """
        Serialize self to a FITS/JSON format.

        Possibility to write it to a new file located at ``path``. Format can be 'fits' or 'json',
        though only the FITS format is officially supported by the IVOA.

        Parameters
        ----------
        path : str, optional
            The path where to save the serialization of self. The serialization is written to ``path`` only if
            ``write_to_file`` is set to True. None by default.
        format : str, optional
            The format of the serialization. Must have its value in ('fits', 'json'). By default, ``format`` is set to
             'fits'.
        optional_kw_dict : dict, optional
            Optional keyword arguments added to the FITS header. Only used if ``format`` is set to 'fits'.
        write_to_file : bool, optional
            Set to False by default. In this case, this method does not write to a file but returns the serialized form
            of self.

        Returns
        -------
        result : `astropy.io.fits.HDUList`/dict
            Depending on the value of ``format``.
        """
        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        # Get all the uniq number from the nuniq intervals
        intervals_uniq_l = IntervalSet.to_nuniq_interval_set(self._interval_set)._intervals
        uniq_l = []
        for uniq_iv in intervals_uniq_l:
            for uniq in range(uniq_iv[0], uniq_iv[1]):
                uniq_l.append(uniq)

        uniq_arr = np.asarray(uniq_l, dtype=np.int64)

        if format == 'fits':
            result = self.__class__._to_fits(uniq_arr=uniq_arr,
                                             moc_order=self.max_order,
                                             optional_kw_dict=optional_kw_dict)
            if write_to_file:
                result.writeto(path, overwrite=True)
        else:
            # json format serialization
            result = self.__class__._to_json(uniq_arr=uniq_arr)
            if write_to_file:
                import json
                with open(path, 'w') as h:
                    h.write(json.dumps(result, sort_keys=True, indent=2))

        return result

    @staticmethod
    def _to_json(uniq_arr):
        """
        Serialize a NUNIQ numpy array to json.

        Parameters
        ----------
        uniq_arr : `~numpy.ndarray`
            An array of uniq numbers corresponding to self.
        Returns
        -------
        result : dict
            A dictionary of HEALPix cells list each indexed by its order.
        """
        result_json = {}

        order_arr, ipix_arr = uniq2orderipix(uniq_arr)
        min_order = order_arr[0]
        max_order = order_arr[-1]

        for order in range(min_order, max_order+1):
            pix_index = np.where(order_arr == order)[0]
            if pix_index.size:
                # there are pixels belonging to the current order
                ipix_order_arr = ipix_arr[pix_index]
                result_json[str(order)] = ipix_order_arr.tolist()

        return result_json

    @staticmethod
    def _to_fits(uniq_arr, moc_order, optional_kw_dict=None):
        """
        Serialize a NUNIQ numpy array to a FITS file.

        Parameters
        ----------
        uniq_arr : `numpy.ndarray`
            An array of uniq numbers corresponding to self.
        moc_order : int
            The maximum order of self.
        optional_kw_dict : dict
            Optional keyword arguments added to the FITS header.

        Returns
        -------
        thdulist : `astropy.io.fits.HDUList`
            The FITS serialization of self.
        """
        if moc_order <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'

        tbhdu = fits.BinTableHDU.from_columns(
            fits.ColDefs([fits.Column(name='UNIQ', format=fits_format, array=uniq_arr)]))
        tbhdu.header['PIXTYPE'] = 'HEALPIX'
        tbhdu.header['ORDERING'] = 'NUNIQ'
        tbhdu.header['COORDSYS'] = ('C', 'reference frame (C=ICRS)')
        tbhdu.header['MOCORDER'] = moc_order
        tbhdu.header['MOCTOOL'] = 'astropy.regions'
        if optional_kw_dict:
            for key in optional_kw_dict:
                tbhdu.header[key] = optional_kw_dict[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def intersection(self, another_moc, *args):
        """
        Intersection between self and other `~regions.MOCSkyRegion` objects.

        Parameters
        ----------
        another_moc : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs

        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object resulting from the intersection of self with the MOCs passed to the method.
        """
        interval_set = self._interval_set.intersection(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.intersection(moc._interval_set)

        return self.__class__(interval_set)

    def union(self, another_moc, *args):
        """
        Union between self and other `~regions.MOCSkyRegion` objects.

        Parameters
        ----------
        another_moc : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs

        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object resulting from the union of self with the MOCs passed to the method.
        """
        interval_set = self._interval_set.union(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.union(moc._interval_set)

        return self.__class__(interval_set)

    def difference(self, another_moc, *args):
        """
        Difference between self and other MOCs.

        Parameters
        ----------
        another_moc : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs
        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object resulting from the difference of self with the MOCs passed to the method.
        """
        interval_set = self._interval_set.difference(another_moc._interval_set)
        for moc in args:
            interval_set = interval_set.difference(moc._interval_set)

        return self.__class__(interval_set)

    def complement(self):
        """
        Compute the complement of self.

        Returns
        -------
        complement : `~regions.MOCSkyRegion`
            A new `~regions.MOCSkyRegion` object corresponding to the complement of self.
        """
        res = []
        intervals_l = sorted(self._interval_set._intervals.tolist())

        if intervals_l[0][0] > 0:
            res.append((0, intervals_l[0][0]))

        last = intervals_l[0][1]

        for itv in intervals_l[1:]:
            res.append((last, itv[0]))
            last = itv[1]

        max_pix_order = 3 << 60

        if last < max_pix_order:
            res.append((last, max_pix_order))

        return self.__class__(IntervalSet(np.asarray(res, dtype=np.int64)))

    def add_neighbours(self):
        """
        Add all the pixels at max order in the neighbourhood of self.

        The algorithm for adding HEALPix cell neighbors follows the steps:

        1. Get the HEALPix cell array of self at its max order.
        2. Get the HEALPix cell array containing the neighbors of the first array (i.e. an ``extended`` HEALPix
           cell array containing the first one).
        3. Subtract the first from the second HEALPix cell array to get only the HEALPix cells
           located at the border of self.
        4. This last HEALPix cell array is added to self.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            Self which has been augmented.
        """
        pix_id_arr = self._best_res_pixels()

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        neighbour_pix_arr = MOCSkyRegion._get_neighbour_pix(hp, pix_id_arr)

        augmented_pix_arr = np.setdiff1d(neighbour_pix_arr, pix_id_arr)

        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - self.max_order)
        intervals_arr = np.vstack((augmented_pix_arr << shift, (augmented_pix_arr + 1) << shift)).T

        self._interval_set = self._interval_set.union(IntervalSet.from_numpy_array(intervals_arr))
        return self

    def remove_neighbours(self):
        """
        Remove all the pixels at max order located at the bound of self.

        The algorithm for removing the HEALPix cells located at the border of self follows the steps:

        1. Get the HEALPix cell array of self at its max order.
        2. Get the HEALPix cell array containing the neighbors of the first array (i.e. an ``extended`` HEALPix
           cell array containing the first one).
        3. Subtract the first from the second HEALPix cell array to get only the HEALPix cells
           located at the border of self.
        4. Same as step 2. to get the HEALPix cell array containing the neighbors of the last computed array (i.e. we get the HEALPix cell neighbors
           of the HEALPix neighbors of self).
        5. This last HEALPix cell array is subtracted from the HEALPix cell array describing self.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            Self whose HEALPix cells located at its border have been removed.
        """
        pix_id_arr = self._best_res_pixels()

        hp = HEALPix(nside=(1 << self.max_order), order='nested')
        neighbour_pix_arr = MOCSkyRegion._get_neighbour_pix(hp, pix_id_arr)

        only_neighbour_arr = np.setxor1d(neighbour_pix_arr, pix_id_arr)

        bound_pix_arr = MOCSkyRegion._get_neighbour_pix(hp, only_neighbour_arr)

        diminished_pix_arr = np.setdiff1d(pix_id_arr, bound_pix_arr)

        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - self.max_order)
        intervals_arr = np.vstack((diminished_pix_arr << shift, (diminished_pix_arr + 1) << shift)).T
        self._interval_set = IntervalSet.from_numpy_array(intervals_arr)
        return self

    @staticmethod
    def _get_neighbour_pix(hp, pix_arr):
        """
        Get the HEALPix cells located in the neighborhood of ``pix_arr``.

        Parameters
        ----------
        hp : `~astropy_healpix.HEALPix`
            The astropy-healpix context.
        pix_arr : `~numpy.ndarray`
            Array of HEALPix cells at a specific order -- linked to ``hp``.

        Returns
        -------
        neighbors_pix_arr : `~numpy.ndarray`
            The HEALPix cells located in the neighborhood of ``pix_arr``.
        """
        neighbors_pix_arr = np.unique(hp.neighbours(pix_arr).ravel())
        # Remove negative HEALPix cell values returned by `~astropy_healpix.HEALPix.neighbours`
        return neighbors_pix_arr[np.where(neighbors_pix_arr >= 0)]

    def _best_res_pixels(self):
        """
        Get the HEALPix cells of self at its maximum order.

        Returns
        -------
        array : `~numpy.ndarray`
            The HEALPix cells of self at its maximum order.
        """
        factor = 2 * (MOCSkyRegion.HPY_MAX_NORDER - self.max_order)
        pix_l = []
        for iv in self._interval_set._intervals:
            for val in range(iv[0] >> factor, iv[1] >> factor):
                pix_l.append(val)

        return np.asarray(pix_l, dtype=np.int64)

    def degrade_to_order(self, new_order):
        """
        Degrade a `~regions.MOCSkyRegion` object.

        The degraded MOC has an order equals to ``new_order``. ``new_order`` must be smaller
        than the current order of the MOC.

        Parameters
        ----------
        new_order : int
            The new maximum order for the degraded MOC

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            The degraded MOC
        """
        shift = 2 * (MOCSkyRegion.HPY_MAX_NORDER - new_order)
        ofs = (int(1) << shift) - 1
        mask = ~ofs
        adda = int(0)
        addb = ofs
        iv_set = []

        for iv in self._interval_set._intervals:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.append((a, b))

        return self.__class__(IntervalSet.from_numpy_array(np.asarray(iv_set, dtype=np.int64)))
