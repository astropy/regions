# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
from .py23_compat import range, int

import numpy as np
from astropy import units as u

from astropy import wcs
from astropy.io import fits
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import Table
from astropy.wcs.utils import pixel_to_skycoord, \
                              skycoord_to_pixel

from astropy_healpix import HEALPix, \
                            lonlat_to_healpix, \
                            nside_to_npix, \
                            level_to_nside, \
                            uniq_to_level_ipix

from .interval_set import IntervalSet
from .utils import trailing_zeros

from ...core import PixCoord, \
                    PixelRegion, \
                    SkyRegion, \
                    BoundingBox, \
                    RegionMask
from ...core.attributes import RegionMeta, \
                               RegionVisual

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
    A MOC (Multi-Order Coverage) map.

    Please refer to the doc example page explaining the basic features
    of this class: :ref:`moc`.

    Parameters
    ----------
    itv_s : `~regions.IntervalSet` object, optional
        A N rows by 2 columns `~numpy.ndarray` storing the set of intervals
        representing the ranges of HEALPix cell indexes of the MOC.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    HPY_MAX_LVL = 29

    def __init__(self, itv_s=None, meta=None, visual=None):
        itv_s = IntervalSet() if itv_s is None else itv_s
        self._itv_s = itv_s
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

    def __eq__(self, other):
        """
        Test equality between two `~regions.MOCSkyRegion`

        Parameters
        ----------
        other : `~regions.MOCSkyRegion`
            The `~regions.MOCSkyRegion` object to test the equality with.

        Returns
        -------
        result : bool
            True is the two `~regions.MOCSkyRegion` are equal.
        """
        if not isinstance(other, MOCSkyRegion):
            raise TypeError('The object you want to test the equality with is not a MOC but a {0}'
                            .format(type(other)))

        return self._itv_s == other._itv_s

    @classmethod
    def from_json(cls, data, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a python dictionary

        Parameters
        ----------
        data : dict
            A dictionary of (str(level), [ipix]) key-value items.
            In the following example, we refer to the MOC containing the [2, 3] HEALPix cells of the level '0'
            plus the [6, 9, 22] cells of the level '5': ``data = {'0': [2, 3], '5': [6, 9, 22]}``
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object.
        """
        itvs = np.array([], dtype=np.int64)

        for level, ipix_l in data.items():
            ipix = np.array(ipix_l, dtype=np.int64)
            p1 = ipix
            p2 = ipix + 1
            shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - int(level))

            itv_batch = np.vstack((p1 << shift, p2 << shift)).T
            if itvs.size == 0:
                itvs = itv_batch
            else:
                itvs = np.vstack((itvs, itv_batch))

        itv_set = IntervalSet(itvs)
        return cls(itv_set, meta, visual)

    @classmethod
    def from_fits(cls, filename, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a FITS file

        Parameters
        ----------
        filename : str
            The path to the FITS file.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object.
        """
        table = Table.read(filename)

        itvs = np.vstack((table['UNIQ'], table['UNIQ']+1)).T

        uniq_itv_s = IntervalSet(itvs)
        itv_s = IntervalSet.from_uniq_itv_s(uniq_itv_s)
        return cls(itv_s, meta, visual)

    @classmethod
    def from_image(cls, header, max_level, mask=None, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from an image stored as a FITS file

        Parameters
        ----------
        header : `~astropy.io.fits.Header`
            The FITS header of the image.
        max_level : int
            The level indicating the maximum resolution of the MOC. Must be <= 29.
        mask : `~numpy.ndarray`, optional
            A 2D boolean numpy array of the same size of the image where pixels evaluated to True are part of
            the new MOC and pixels evaluated to False are not.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object.
        """
        # Get the image dimensions
        height = header['NAXIS2']
        width = header['NAXIS1']

        # Get the image WCS to retrieve back its world coordinates
        w = wcs.WCS(header)

        # Get the pixels coordinates
        if mask is not None:
            # Get the pixels in the mask if available
            y, x = np.where(mask)
            pix = np.dstack((x, y))[0]
        else:
            # Get a uniform grid of pixels whose boundaries are the dimensions of the image.
            step = 1
            x, y = np.mgrid[0.5:(width + 0.5 + step):step, 0.5:(height + 0.5 + step):step]
            pix = np.dstack((x.ravel(), y.ravel()))[0]

        # Retrieve back the world coordinates using the image WCS
        world_pix = w.wcs_pix2world(pix, 1)

        # Call the from_lonlat method to create the MOC from the world coordinates.
        return MOCSkyRegion.from_lonlat(lon=world_pix[:, 0] * u.deg,
         lat=world_pix[:, 1] * u.deg,
         max_level=max_level,
         meta=meta,
         visual=visual)

    @classmethod
    def from_skycoord(cls, skycoord, max_level, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from a `~astropy.coordinates.SkyCoord` object

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            An astropy `~astropy.coordinates.SkyCoord` object.
        max_level : int
            The level indicating the maximum resolution of the MOC. Must be <= 29.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        nside = level_to_nside(max_level)
        hp = HEALPix(nside=nside, order='nested')
        ipix = hp.lonlat_to_healpix(skycoord.icrs.ra, skycoord.icrs.dec)

        shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - max_level)
        itvs = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return cls(itv_s, meta, visual)

    @classmethod
    def from_lonlat(cls, lon, lat, max_level, meta=None, visual=None):
        """
        Create a `~regions.MOCSkyRegion` from ``lon``, ``lat`` `astropy.units.Quantity`.

        Parameters
        ----------
        lon : `astropy.units.Quantity`
            A set of astropy formatted ra quantities.
        lat : `astropy.units.Quantity`
            A set of astropy formatted dec quantities.
        max_level : int
            The level indicating the maximum resolution of the MOC. Must be <= 29.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        nside = level_to_nside(max_level)
        hp = HEALPix(nside=nside, order='nested')
        ipix = hp.lonlat_to_healpix(lon, lat)

        shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - max_level)
        itvs = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return cls(itv_s, meta, visual)

    @property
    def max_level(self):
        """
        Returns the level of the smallest HEALPix cells
        """
        combo = int(0)
        for iv in self._itv_s._data:
            combo |= iv[0] | iv[1]

        ret = MOCSkyRegion.HPY_MAX_LVL - (trailing_zeros(combo) // 2)
        if ret < 0:
            ret = 0

        return ret

    @property
    def sky_fraction(self):
        """
        Returns the sky fraction coverage percentage
        """
        # Get the number of ipixels at the deepest level
        ipix = self._best_res_pixels()
        num_ipix = ipix.size

        # Get the maximum number of ipixels for the level of the deepest ipixels
        npix = 3 << (2*(self.max_level + 1))

        return num_ipix / float(npix)

    def contains(self, ra, dec):
        """
        Check whether some sky positions fall inside the MOC

        The size of ``ra`` and ``dec`` must be equal.

        Parameters
        ----------
        ra : `~astropy.units.Quantity`
            A set of astropy formatted ra quantities
        dec: `~astropy.units.Quantity`
            A set of astropy formatted dec quantities

        Returns
        -------
        array : `~numpy.ndarray`
            A boolean numpy array telling which positions are inside the MOC depending on the
            value of ``self.meta['include']``.
        """
        # Build the mask boolean array of the ipixels from the MOC
        # at the deepest level.
        nside = level_to_nside(self.max_level)
        npix = nside_to_npix(nside)
        mask = np.zeros(npix, dtype=bool)

        ipix = self._best_res_pixels()
        mask[ipix] = True

        if self.meta.get('include', False):
            mask = np.logical_not(mask)

        # Retrieve the ipixels containing the sky positions
        hp = HEALPix(nside=nside, order='nested')
        ipix = hp.lonlat_to_healpix(ra, dec)

        return mask[ipix]

    def write(self, path, format='fits', fits_optional_kw=None):
        """
        Write to a FITS/JSON file.

        Parameters
        ----------
        path : str
            The path where to save the file.
        format : str, optional
            The format of the serialization. Must have its value in ('fits', 'json'). By default, ``format`` is set to
             'fits' as FITS is the one officially supported by the IVOA.
        fits_optional_kw : dict, optional
            Optional arguments added to the FITS header. Only used if ``format`` is set to 'fits'.
        """
        serialization = self.serialize(format=format, fits_optional_kw=fits_optional_kw)

        if format == 'fits':
            serialization.writeto(path, overwrite=True)
        else:
            # json format serialization
            import json
            with open(path, 'w') as h:
                h.write(json.dumps(serialization, sort_keys=True, indent=2))

    def serialize(self, format='fits', fits_optional_kw=None):
        """
        Serialize to a FITS/JSON format

        Parameters
        ----------
        format : str, optional
            The format of the serialization. Must have its value in ('fits', 'json'). By default, ``format`` is set to
             'fits' as FITS is the one officially supported by the IVOA.
        fits_optional_kw : dict, optional
            Optional arguments added to the FITS header. Only used if ``format`` is set to 'fits'.

        Returns
        -------
        result : `astropy.io.fits.HDUList`/dict
            Depending whether the specified format is FITS (resp. JSON).
        """
        fmts = ('fits', 'json')
        if format not in fmts:
            raise ValueError('format should be one of %s' % (str(fmts)))

        # Get all the uniq number from the nuniq intervals
        uniq_itv_s = IntervalSet.to_uniq_itv_s(self._itv_s)
        uniq_itvs = uniq_itv_s._data
        uniq_pix_l = []
        for uniq_itv in uniq_itvs:
            for uniq in range(uniq_itv[0], uniq_itv[1]):
                uniq_pix_l.append(uniq)

        uniq_pix = np.asarray(uniq_pix_l, dtype=np.int64)

        if format == 'fits':
            result = self._to_fits(uniq_pix=uniq_pix,
                                   fits_optional_kw=fits_optional_kw)
        else:
            # JSON format serialization
            result = self._to_json(uniq_pix=uniq_pix)

        return result

    def _to_json(self, uniq_pix):
        """
        Serialize a numpy array of uniq numbers to JSON

        Parameters
        ----------
        uniq_pix : `~numpy.ndarray`
            A numpy array of uniq numbers.
        Returns
        -------
        res : dict
            A python dict containing (level, [ipix]) key-value items
        """
        res = {}

        levels, ipix = uniq_to_level_ipix(uniq_pix)
        min_lvl = levels[0]
        max_lvl = levels[-1]

        for lvl in range(min_lvl, max_lvl+1):
            index = np.where(levels == lvl)[0]
            if ipix.size:
                # There are ipixels of this specific level
                ipix_lvl = ipix[index]
                res[str(lvl)] = ipix_lvl.tolist()

        return res

    def _to_fits(self, uniq_pix, fits_optional_kw=None):
        """
        Serialize a numpy array of uniq numbers to FITS

        Parameters
        ----------
        uniq_pix : `numpy.ndarray`
            A numpy array of uniq numbers.
        fits_optional_kw : dict
            Optional arguments added to the FITS header.

        Returns
        -------
        thdulist : `astropy.io.fits.HDUList`
            The FITS serialization of self
        """
        max_level = self.max_level
        if max_level <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'

        tbhdu = fits.BinTableHDU.from_columns(
            fits.ColDefs([fits.Column(name='UNIQ', format=fits_format, array=uniq_pix)]))
        tbhdu.header['PIXTYPE'] = 'HEALPIX'
        tbhdu.header['ORDERING'] = 'NUNIQ'
        tbhdu.header['COORDSYS'] = ('C', 'reference frame (C=ICRS)')
        tbhdu.header['MOCORDER'] = max_level
        tbhdu.header['MOCTOOL'] = 'astropy.regions'
        if fits_optional_kw:
            for key in fits_optional_kw:
                tbhdu.header[key] = fits_optional_kw[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def intersection(self, other, *args):
        """
        Intersection with other `~regions.MOCSkyRegion` objects

        Parameters
        ----------
        other : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs

        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        itv = self._itv_s.intersection(other._itv_s)
        for moc in args:
            itv = itv.intersection(moc._itv_s)

        return self.__class__(itv)

    def union(self, other, *args):
        """
        Union with other `~regions.MOCSkyRegion` objects

        Parameters
        ----------
        other : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs

        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        itv = self._itv_s.union(other._itv_s)
        for moc in args:
            itv = itv.union(moc._itv_s)

        return self.__class__(itv)

    def difference(self, other, *args):
        """
        Difference with other `~regions.MOCSkyRegion` objects

        Parameters
        ----------
        other : `~regions.MOCSkyRegion`
            Another mandatory MOC
        args : `~regions.MOCSkyRegion`, optional
            More MOCs
        Returns
        -------
        result : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        itv = self._itv_s.difference(other._itv_s)
        for moc in args:
            itv = itv.difference(moc._itv_s)

        return self.__class__(itv)

    def complement(self):
        """
        Return the MOC complement

        Returns
        -------
        complement : `~regions.MOCSkyRegion`
            The complement given as a `~regions.MOCSkyRegion` object
        """
        res = []
        itvs_l = sorted(self._itv_s._data.tolist())

        if itvs_l[0][0] > 0:
            res.append((0, itvs_l[0][0]))

        last = itvs_l[0][1]

        for itv in itvs_l[1:]:
            res.append((last, itv[0]))
            last = itv[1]

        max_npix = 3 << 60

        if last < max_npix:
            res.append((last, max_npix))

        itvs = np.asarray(res, dtype=np.int64)
        itv_s = IntervalSet(itvs)
        return self.__class__(itv_s)

    def add_neighbours(self):
        """
        Add the HEALPix cell indexes located at the border

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            The augmented MOC
        """
        # Retrieve the ipixels at the deepest level
        ipix = self._best_res_pixels()
        max_level = self.max_level
        nside = level_to_nside(max_level)

        # Compute the neighbours of the ipixels retrieved
        hp = HEALPix(nside=nside, order='nested')
        ipix_neigh = hp.neighbours(ipix)[[0, 2, 4, 6], :]

        # Get the union of the ipixels with their neighbours
        res = np.union1d(ipix, ipix_neigh)

        shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - max_level)
        itvs = np.vstack((res << shift, (res + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return MOCSkyRegion(itv_s)

    def remove_neighbours(self):
        """
        Remove the HEALPix cell indexes located at the border

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            The diminished MOC
        """
        # Retrieve the ipixels at the deepest level
        ipix = self._best_res_pixels()
        max_level = self.max_level
        nside = level_to_nside(max_level)

        # Retrieve the ipixels being at the border
        hp = HEALPix(nside=nside, order='nested')
        ipix_neigh = hp.neighbours(ipix)[[0, 2, 4, 6], :]
        mask = np.isin(ipix_neigh, ipix)
        num_neigh = mask.sum(axis=0)
        border = num_neigh < 4

        # Get the ipixels which are not at the border
        res = ipix[~border]

        shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - max_level)
        itvs = np.vstack((res << shift, (res + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return MOCSkyRegion(itv_s)

    def _best_res_pixels(self):
        """
        Return the HEALPix cell indexes at its deepest level.

        Returns
        -------
        array : `~numpy.ndarray`
            A numpy array containing the HEALPix cell indexes at the deepest level.
        """
        factor = 2 * (MOCSkyRegion.HPY_MAX_LVL - self.max_level)
        ipix = []
        for iv in self._itv_s._data:
            for val in range(iv[0] >> factor, iv[1] >> factor):
                ipix.append(val)

        return np.asarray(ipix, dtype=np.int64)

    def degrade_to_order(self, new_level):
        """
        Degrade to a shallower level

        The degraded MOC will have a maximum level equal to ``new_level``. ``new_level`` must be smaller
        than the current max level of the MOC.

        Parameters
        ----------
        new_level : int
            The new maximum level for the degraded MOC

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            A `~regions.MOCSkyRegion` object
        """
        shift = 2 * (MOCSkyRegion.HPY_MAX_LVL - new_level)
        ofs = (int(1) << shift) - 1
        mask = ~ofs
        adda = int(0)
        addb = ofs
        itvs = []

        for iv in self._itv_s._data:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                itvs.append((a, b))

        itvs = np.asarray(itvs, dtype=np.int64)
        itv_s = IntervalSet(itvs)
        return self.__class__(itv_s)
