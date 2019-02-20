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
                            nside_to_npix, \
                            level_to_nside, \
                            uniq_to_level_ipix

from .interval_set import IntervalSet
from .utils import trailing_zeros

from .boundaries import Boundaries
from .plot import fill
from .plot import border
from .plot import axis_viewport
from .plot import culling_backfacing_cells

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
    A MOC (Multi-Order Coverage) map in pixel coordinates.

    Parameters
    ----------
    sky_region : `~regions.MOCSkyRegion`
        The `~regions.MOCSkyRegion` instance.
    meta : `~regions.RegionMeta`, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    def __init__(self, wcs, sky_region, meta=None, visual=None):
        self.moc = sky_region
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('_interval_set',)
        self.wcs = wcs
        self.vertices = self._compute_vertices()

    def _compute_vertices(self):
        # Simplify the MOC for plotting purposes:
        # 1. Degrade the MOC if the FOV is enough big so that we cannot see the smallest HEALPix cells.
        # 2. For small FOVs, plot the MOC & POLYGONAL_MOC_FROM_FOV
        moc = fill.build_plotting_moc(moc=self.moc, wcs=self.wcs)
        # If the FOV contains no cells, then moc_to_plot (i.e. the intersection between the moc
        # and the MOC created from the FOV polygon) will be empty.
        # If it is the case, we exit the method without doing anything.
        vertices = np.array([], dtype=float)
        if moc.empty():
            return vertices

        depth_ipix_d = moc.serialize(format="json")
        depth_ipix_clean_d = culling_backfacing_cells.from_moc(depth_ipix_d=depth_ipix_d, wcs=self.wcs)

        for depth, ipix in depth_ipix_clean_d.items():
            step = 1
            depth = int(depth)
            if depth < 3:
                step = 2

            nside = level_to_nside(depth)
            hp = HEALPix(order="nested", nside=nside, frame=ICRS())
            ipix_boundaries = hp.boundaries_skycoord(ipix, step=step)
            # Projection on the given WCS
            xp, yp = skycoord_to_pixel(ipix_boundaries, wcs=self.wcs)

            patches = np.vstack((xp.flatten(), yp.flatten())).T
            patches = patches.reshape((-1, step*4, 2))

            if vertices.size == 0:
                vertices = patches
            else:
                vertices = np.append(vertices, patches, axis=0)

        return vertices

    def as_artist(self, origin=(0, 0), **kw_mpl_pathpatch):
        """
        Matplotlib patch object for this region (`matplotlib.patches.PathPatch`)

        Parameters:
        -----------
        origin : `~numpy.ndarray`, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs: dict
            All keywords that a `~matplotlib.patches.PathPatch` object accepts

        Returns
        -------
        patch : `~matplotlib.patches.PathPatch`
            Matplotlib path patch
        """
        pathpatch_mpl = fill.fill(moc=self.moc, wcs=self.wcs, origin=origin, **kw_mpl_pathpatch)
        return pathpatch_mpl

    def plot(self, ax=None, origin=(0, 0), **kwargs):
        """
        Calls ``as_artist`` method forwarding all kwargs and adds patch
        to given axis.

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        ax : `~matplotlib.axes.Axes`, optional
            Axis
        perimeter : bool, optional
            Plot the perimeter of the MOCPixelRegion. Default to True.
        kwargs : `dict`
            keywords that a `~matplotlib.patches.Patch` accepts

        Returns
        -------
        ax : `~matplotlib.axes.Axes`
            Axes on which the patch is added.
        """
        PixelRegion.plot(self, ax=ax, origin=origin, **kwargs)
        axis_viewport.set(ax=ax, wcs=self.wcs)
        return ax

    def to_sky(self, wcs):
        return self.moc

    @property
    def bounding_box(self):
        xmin = np.min(self.vertices[:, :, 0])
        xmax = np.max(self.vertices[:, :, 0])
        ymin = np.min(self.vertices[:, :, 1])
        ymax = np.max(self.vertices[:, :, 1])
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

        vx, vy = self.vertices[:, :, 0], self.vertices[:, :, 1]

        # Loop over all the projeted HEALPix cells to get their overlap grids
        buffer = np.zeros(shape=(ny, nx))
        for i in range(vx.shape[0]):
            ovlp_grid = polygonal_overlap_grid(
                xmin, xmax, ymin, ymax,
                nx, ny, vx[i], vy[i],
                use_exact, subpixels,
            )
            # Sum this overlap grid to the buffer
            buffer += ovlp_grid

        # Clip its values to the interval (0, 1)
        mask = np.clip(a=buffer, a_min=0, a_max=1)
        return RegionMask(mask, bbox=bbox)

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
        vx, vy = self.vertices[:, :, 0], self.vertices[:, :, 1]
        mask = points_in_polygons(x.flatten(), y.flatten(), vx, vy).astype(bool)
        in_poly = mask.reshape(shape)
        if self.meta.get('include', True):
            return in_poly
        else:
            return np.logical_not(in_poly)

class MOCSkyRegion(SkyRegion):
    """
    Multi-order spatial coverage class

    A MOC describes the coverage of an arbitrary region on the unit sphere.
    MOCs are usually used for describing the global coverage of catalog/image surveys such as GALEX or SDSS.
    A MOC corresponds to a list of `HEALPix <https://healpix.sourceforge.io/>`__ cells at different depths.
    This class gives you the possibility to:

    1. Define `regions.MOCSkyRegion` objects:

    - From a FITS file that stores HEALPix cells (see `from_fits`).
    - Directly from a list of HEALPix cells expressed either as a numpy structural array (see `from_cells`) or a simple
      python dictionnary (see `from_json`).
    - From a list of sky coordinates (see `from_skycoord`, `from_lonlat`).
    - From a convex/concave polygon (see `from_polygon`).
    - From a cone (will be implemented in a next version).

    2. Perform fast logical operations between `regions.MOCSkyRegion` objects:

    - The `intersection`
    - The `union`
    - The `difference`
    - The `complement`

    3. Plot the `regions.MOCSkyRegion` objects:

    - Draw the MOC with its HEALPix cells (see `fill`)
    - Draw the perimeter of a MOC (see `border`)

    4. Get the sky coordinates defining the border(s) of `regions.MOCSkyRegion` objects (see `get_boundaries`).

    5. Serialize `regions.MOCSkyRegion` objects to `astropy.io.fits.HDUList` or JSON dictionary and save it to a file.

    Parameters
    ----------
    itv_s : `~regions.IntervalSet` object, optional
        A N rows by 2 columns `~numpy.ndarray` storing the ranges of HEALPix cells contained in the MOC.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    HPY_MAX_DEPTH = 29

    def __init__(self, itv=None, meta=None, visual=None):
        interval = IntervalSet() if itv is None else itv
        self._itv = interval
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        """
        Convert the `~regions.MOCSkyRegion` instance to a `~regions.MOCPixelRegion`

        The MOC is projected into the image coordinates system using an `~astropy.wcs.WCS` instance.
        The HEALPix cells that backface the projection are removed.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The world coordinate system <-> Image coordinate system projection

        Returns
        -------
        pixel_region : `~regions.MOCPixelRegion`
            The projection of the MOCSkyRegion instance.
        """
        return MOCPixelRegion(wcs=wcs, sky_region=self, meta=self.meta, visual=self.visual)

    def __repr__(self):
        return self._itv.__repr__()

    def __eq__(self, other):
        """
        Test equality between two `~regions.MOCSkyRegion`

        Parameters
        ----------
        other : `~regions.MOCSkyRegion`
            The `~regions.MOCSkyRegion` instance to test the equality with.

        Returns
        -------
        result : bool
            True is the two `~regions.MOCSkyRegion` are equal.
        """
        if not isinstance(other, MOCSkyRegion):
            raise TypeError('Cannot compare a MOCSkyRegion with a {0}'.format(type(other)))

        return self._itv == other._itv

    def empty(self):
        """
        Checks whether the `~regions.MOCSkyRegion` instance is empty

        A MOC is empty when it contains no ranges of HEALPix cell indexes.

        Returns
        -------
        result: bool
            True if the MOC instance is empty.
        """
        return self._itv.empty()

    @property
    def max_depth(self):
        """
        Depth of the smallest HEALPix cells found in the MOC instance.
        """
        combo = int(0)
        for iv in self._itv._data:
            combo |= iv[0] | iv[1]

        ret = MOCSkyRegion.HPY_MAX_DEPTH - trailing_zeros(combo)//2
        if ret < 0:
            ret = 0

        return ret

    def intersection(self, other, *args):
        """
        Intersection between the MOC instance and other MOCs.

        Parameters
        ----------
        other : `regions.MOCSkyRegion`
            The MOC used for performing the intersection with self.
        args : `regions.MOCSkyRegion`
            Other additional MOCs to perform the intersection with.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC.
        """
        itv = self._itv.intersection(other._itv)
        for moc in args:
            itv = itv.intersection(moc._itv)

        return self.__class__(itv)

    def union(self, other, *args):
        """
        Union between the MOC instance and other MOCs.

        Parameters
        ----------
        other : `regions.MOCSkyRegion`
            The MOC used for performing the union with self.
        args : `regions.MOCSkyRegion`
            Other additional MOCs to perform the union with.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC.
        """
        itv = self._itv.union(other._itv)
        for moc in args:
            itv = itv.union(moc._itv)

        return self.__class__(itv)

    def difference(self, other, *args):
        """
        Difference between the MOC instance and other MOCs.

        Parameters
        ----------
        other : `regions.MOCSkyRegion`
            The MOC used that will be substracted to self.
        args : `regions.MOCSkyRegion`
            Other additional MOCs to perform the difference with.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC.
        """
        itv = self._itv.difference(other._itv)
        for moc in args:
            itv = itv.difference(moc._itv)

        return self.__class__(itv)

    def complement(self):
        """
        Returns the complement of the MOC instance.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC.
        """
        res = []
        itvs_l = sorted(self._itv._data.tolist())

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
        itv = IntervalSet(itvs)
        return self.__class__(itv)


    @classmethod
    def from_cells(cls, cells, meta=None, visual=None):
        """
        Creates a MOC from a numpy array representing the HEALPix cells.

        Parameters
        ----------
        cells : `numpy.ndarray`
            Must be a numpy structured array (See https://docs.scipy.org/doc/numpy-1.15.0/user/basics.rec.html).
            The structure of a cell contains 3 attributes:

            - A `ipix` value being a np.uint64
            - A `depth` value being a np.uint32
            - A `fully_covered` flag bit stored in a np.uint8
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            The MOC.
        """
        shift = (MOCSkyRegion.HPY_MAX_DEPTH - cells["depth"]) << 1

        p1 = cells["ipix"]
        p2 = cells["ipix"] + 1

        itvs = np.vstack((p1 << shift, p2 << shift)).T

        return cls(IntervalSet(itvs), meta, visual)

    @classmethod
    def from_json(cls, data, meta=None, visual=None):
        """
        Creates a `~regions.MOCSkyRegion` instance from a dictionary of HEALPix cell arrays indexed by their depth.

        Parameters
        ----------
        data : dict(str : [int])
            A dictionary of HEALPix cell arrays indexed by their depth.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            The MOC instance.

        Examples
        --------
        >>> from regions import MOCSkyRegion
        >>> moc = MOCSkyRegion.from_json({
        ...     "0": [0, 3, 4],
        ...     "1": [5, 6],
        ... })
        """
        itvs = np.array([], dtype=np.int64)
        for order, pix_l in data.items():
            if len(pix_l) == 0:
                continue
            pix = np.array(pix_l, dtype=np.int64)
            p1 = pix
            p2 = pix + 1
            shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - int(order))

            itv = np.vstack((p1 << shift, p2 << shift)).T
            if itvs.size == 0:
                itvs = itv
            else:
                itvs = np.vstack((itvs, itv))

        itv_s = IntervalSet(itvs)
        return cls(itv_s, meta, visual)

    def _uniq_pixels_iterator(self):
        """
        Generator giving the NUNIQ HEALPix pixels of the MOC.

        Returns
        -------
        uniq :
            the NUNIQ HEALPix pixels iterator
        """
        itvs_uniq_l = IntervalSet.to_uniq_itv_s(self._itv)._data
        for uniq_iv in itvs_uniq_l:
            for uniq in np.arange(uniq_iv[0], uniq_iv[1], dtype=np.int64):
                yield uniq

    @classmethod
    def from_fits(cls, filename, meta=None, visual=None):
        """
        Creates a `~regions.MOCSkyRegion` from a FITS file.

        The specified FITS file must store the MOC (i.e. the list of HEALPix cells it contains) in a binary HDU table.

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
            The MOC instance.
        """
        table = Table.read(filename)

        itvs = np.vstack((table['UNIQ'], table['UNIQ']+1)).T

        nuniq_itv = IntervalSet(itvs)
        itv = IntervalSet.from_uniq_itv_s(nuniq_itv)
        return cls(itv, meta, visual)

    @staticmethod
    def _to_json(uniq):
        """
        Serializes a MOC to the JSON format.

        Parameters
        ----------
        uniq : `~numpy.ndarray`
            The array of HEALPix cells representing the MOC to serialize.

        Returns
        -------
        result : dict(str : [int])
            A dictionary of HEALPix cell lists indexed by their depth.
        """
        result = {}

        depth, ipix = uniq_to_level_ipix(uniq)
        min_depth = np.min(depth[0])
        max_depth = np.max(depth[-1])

        for d in range(min_depth, max_depth+1):
            pix_index = np.where(depth == d)[0]
            if pix_index.size:
                # There are pixels belonging to the current order
                ipix_depth = ipix[pix_index]
                result[str(d)] = ipix_depth.tolist()

        return result

    def _to_fits(self, uniq_pix, optional_kw_dict=None):
        """
        Serializes a MOC to the FITS format.

        Parameters
        ----------
        uniq_pix : `numpy.ndarray`
            The array of HEALPix cells representing the MOC to serialize.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header.

        Returns
        -------
        thdulist : `astropy.io.fits.HDUList`
            The list of HDU tables.
        """
        depth = self.max_depth
        if depth <= 13:
            fits_format = '1J'
        else:
            fits_format = '1K'

        tbhdu = fits.BinTableHDU.from_columns(
            fits.ColDefs([
                fits.Column(name='UNIQ', format=fits_format, array=uniq_pix)
            ]))
        tbhdu.header['PIXTYPE'] = 'HEALPIX'
        tbhdu.header['ORDERING'] = 'NUNIQ'
        tbhdu.header['COORDSYS'] = ('C', 'reference frame (C=ICRS)')
        tbhdu.header['MOCORDER'] = depth
        tbhdu.header['MOCTOOL'] = 'MOCPy'
        if optional_kw_dict:
            for key in optional_kw_dict:
                tbhdu.header[key] = optional_kw_dict[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def serialize(self, format='fits', optional_kw_dict=None):
        """
        Serializes the MOC into a specific format.

        Possible formats are FITS and JSON.

        Parameters
        ----------
        format : str
            'fits' by default. The other possible choice is 'json'.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.

        Returns
        -------
        result : `astropy.io.fits.HDUList` or JSON dictionary
            The result of the serialization.
        """
        formats = ('fits', 'json')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        uniq_l = []
        for uniq in self._uniq_pixels_iterator():
            uniq_l.append(uniq)

        uniq = np.array(uniq_l, dtype=np.int64)

        if format == 'fits':
            result = self._to_fits(uniq_pix=uniq,
                                   optional_kw_dict=optional_kw_dict)
        else:
            # json format serialization
            result = self.__class__._to_json(uniq=uniq)

        return result

    def write(self, path, format='fits', optional_kw_dict=None):
        """
        Writes the MOC to a file.

        Format can be 'fits' or 'json', though only the fits format is officially supported by the IVOA.

        Parameters
        ----------
        path : str, optional
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC will be serialized before being saved. Possible formats are "fits" or "json".
            By default, ``format`` is set to "fits".
        optional_kw_dict : optional
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.
        """
        serialization = self.serialize(format=format, optional_kw_dict=optional_kw_dict)
        if format == 'fits':
            serialization.writeto(path, overwrite=True)
        else:
            import json
            with open(path, 'w') as h:
                h.write(json.dumps(serialization, sort_keys=True, indent=2))

    def degrade_to_depth(self, new_depth):
        """
        Degrades the MOC instance to a new, less precise, MOC.

        The maximum depth (i.e. the depth of the smallest HEALPix cells that can be found in the MOC) of the
        degraded MOC is set to ``new_depth``.

        Parameters
        ----------
        new_depth : int

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            The degraded MOC.
        """
        shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - new_depth)
        ofs = (int(1) << shift) - 1
        mask = ~ofs
        adda = int(0)
        addb = ofs
        iv_set = []

        for iv in self._itv._data:
            a = (iv[0] + adda) & mask
            b = (iv[1] + addb) & mask
            if b > a:
                iv_set.append((a, b))

        return self.__class__(IntervalSet(np.asarray(iv_set, dtype=np.int64)))

    def _best_res_pixels(self):
        """
        Returns a numpy array of all the HEALPix indexes contained in the MOC at its max depth.

        Returns
        -------
        result : `~numpy.ndarray`
            The array of HEALPix at ``max_depth``
        """
        factor = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - self.max_depth)
        pix_l = []
        for iv in self._itv._data:
            for val in np.arange(iv[0] >> factor, iv[1] >> factor, dtype=np.int64):
                pix_l.append(val)

        return np.asarray(pix_l, dtype=np.int64)

    def contains(self, ra, dec, keep_inside=True):
        """
        Returns a boolean mask array of the positions lying inside (or outside) the MOC instance.

        Parameters
        ----------
        ra : `astropy.units.Quantity`
            Right ascension array
        dec : `astropy.units.Quantity`
            Declination array

        Returns
        -------
        array : `~numpy.ndarray`
            A boolean numpy array telling which positions are inside the MOC depending on the
            value of ``self.meta['include']``.
        """
        depth = self.max_depth

        nside = level_to_nside(depth)
        npix = nside_to_npix(nside)
        mask = np.zeros(npix, dtype=bool)

        ipix = self._best_res_pixels()
        mask[ipix] = True

        if self.meta.get('include', False):
            mask = np.logical_not(mask)

        # Retrieve the cells containing the sky positions
        hp = HEALPix(nside=nside, order='nested')
        ipix = hp.lonlat_to_healpix(ra, dec)

        return mask[ipix]

    def add_neighbours(self):
        """
        Extends the MOC instance so that it includes the HEALPix cells touching its border.

        The depth of the HEALPix cells added at the border is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            A new MOC instance that extends of one degree the initial MOC.
        """
        # Retrieve the ipixels at the deepest level
        ipix = self._best_res_pixels()
        max_depth = self.max_depth
        nside = level_to_nside(max_depth)

        # Compute the neighbours of the ipixels retrieved
        hp = HEALPix(nside=nside, order='nested')
        ipix_neigh = hp.neighbours(ipix)[[0, 2, 4, 6], :]
        ipix_neigh = ipix_neigh.reshape((1, -1))[0]

        # Get the union of the ipixels with their neighbours
        res = np.union1d(ipix, ipix_neigh)

        shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - max_depth)
        itvs = np.vstack((res << shift, (res + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return MOCSkyRegion(itv_s)

    def remove_neighbours(self):
        """
        Removes from the MOC instance the HEALPix cells located at its border.

        The depth of the HEALPix cells removed is equal to the maximum depth of the MOC instance.

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            A new MOC instance from which the HEALPix cells located at the border
            of the initial MOC have been removed.
        """
        # Retrieve the ipixels at the deepest level
        ipix = self._best_res_pixels()
        max_depth = self.max_depth
        nside = level_to_nside(max_depth)

        # Retrieve the ipixels being at the border
        hp = HEALPix(nside=nside, order='nested')
        ipix_neigh = hp.neighbours(ipix)[[0, 2, 4, 6], :]

        r1 = np.in1d(ipix_neigh[0, :], ipix)
        r2 = np.in1d(ipix_neigh[1, :], ipix)
        r3 = np.in1d(ipix_neigh[2, :], ipix)
        r4 = np.in1d(ipix_neigh[3, :], ipix)

        mask = np.vstack((r1, r2, r3, r4))

        num_neigh = mask.sum(axis=0)
        border = num_neigh < 4

        # Get the ipixels which are not at the border
        res = ipix[~border]

        shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - max_depth)
        itvs = np.vstack((res << shift, (res + 1) << shift)).T

        itv_s = IntervalSet(itvs)
        return MOCSkyRegion(itv_s)

    def fill(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draws the MOC on a matplotlib axis.

        This performs the projection of the cells from the world coordinate system to the pixel image coordinate system.
        You are able to specify various styling kwargs for `matplotlib.patches.PathPatch`
        (see the `list of valid keywords <https://matplotlib.org/api/_as_gen/matplotlib.patches.PathPatch.html#matplotlib.patches.PathPatch>`__).

        Parameters
        ----------
        ax : `matplotlib.axes.Axes`
            Matplotlib axis.
        wcs : `astropy.wcs.WCS`
            WCS defining the World system <-> Image system projection.
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`.

        Examples
        --------
        >>> from astropy.utils.data import get_pkg_data_filename
        >>> from regions import MOCSkyRegion, WCS
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions')
        >>> moc = MOCSkyRegion.from_fits(filename)
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with WCS(fig,
        ...         fov=50 * u.deg,
        ...         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
        ...         coordsys="icrs",
        ...         rotation=Angle(0, u.degree),
        ...         projection="AIT") as wcs:
        ...     ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ...     # Call fill giving the matplotlib axe and the `~astropy.wcs.WCS` object.
        ...     # We will set the matplotlib keyword linewidth to 0 so that it does not plot
        ...     # the border of each HEALPix cell.
        ...     # The color can also be specified along with an alpha value.
        ...     moc.fill(ax=ax, wcs=wcs, linewidth=0, alpha=0.5, fill=True, color="green")
        >>> plt.xlabel('ra')
        >>> plt.ylabel('dec')
        >>> plt.grid(color="black", linestyle="dotted")
        """
        mpl_pathpatch = fill.fill(moc=self, wcs=wcs, origin=(0, 0), **kw_mpl_pathpatch)
        ax.add_patch(mpl_pathpatch)
        axis_viewport.set(ax=ax, wcs=wcs)

    def border(self, ax, wcs, **kw_mpl_pathpatch):
        """
        Draws the MOC border(s) on a matplotlib axis.

        This performs the projection of the sky coordinates defining the perimeter of the MOC to the pixel image coordinate system.
        You are able to specify various styling kwargs for `matplotlib.patches.PathPatch`
        (see the `list of valid keywords <https://matplotlib.org/api/_as_gen/matplotlib.patches.PathPatch.html#matplotlib.patches.PathPatch>`__).

        Parameters
        ----------
        ax : `matplotlib.axes.Axes`
            Matplotlib axis.
        wcs : `astropy.wcs.WCS`
            WCS defining the World system <-> Image system projection.
        kw_mpl_pathpatch
            Plotting arguments for `matplotlib.patches.PathPatch`

        Examples
        --------
        >>> from astropy.utils.data import get_pkg_data_filename
        >>> from regions import MOCSkyRegion, WCS
        >>> from astropy.coordinates import Angle, SkyCoord
        >>> import astropy.units as u
        >>> # Load a MOC, e.g. the MOC of GALEXGR6-AIS-FUV
        >>> filename = get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions')
        >>> moc = MOCSkyRegion.from_fits(filename)
        >>> # Plot the MOC using matplotlib
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure(111, figsize=(15, 15))
        >>> # Define a WCS as a context
        >>> with WCS(fig,
        ...         fov=50 * u.deg,
        ...         center=SkyCoord(0, 20, unit='deg', frame='icrs'),
        ...         coordsys="icrs",
        ...         rotation=Angle(0, u.degree),
        ...         projection="AIT") as wcs:
        ...     ax = fig.add_subplot(1, 1, 1, projection=wcs)
        ...     # Call border giving the matplotlib axe and the `~astropy.wcs.WCS` object.
        ...     moc.border(ax=ax, wcs=wcs, alpha=0.5, color="red")
        >>> plt.xlabel('ra')
        >>> plt.ylabel('dec')
        >>> plt.grid(color="black", linestyle="dotted")
        """
        border.border(self, ax, wcs, **kw_mpl_pathpatch)

    def get_boundaries(self, depth=None):
        """
        Returns the sky coordinates defining the border(s) of the MOC.

        The border(s) are expressed as a list of SkyCoord.
        Each SkyCoord refers to the coordinates of one border of the MOC (i.e.
        either a border of a connexe MOC part or a border of a hole
        located in a connexe MOC part).

        Parameters
        ----------
        depth : int
            The depth of the MOC before computing its boundaries.
            A shallow depth leads to a faster computation.
            By default the maximum depth of the MOC is taken.

        Return
        ------
        coords: [`~astropy.coordinates.SkyCoord`]
            A list of `~astropy.coordinates.SkyCoord` each describing one border.
        """
        return Boundaries.get(self, depth)

    @classmethod
    def from_image(cls, header, max_depth, mask=None, meta=None, visual=None):
        """
        Creates a `~regions.MOCSkyRegion` instance from an image stored as a FITS file.

        Parameters
        ----------
        header : `astropy.io.fits.Header`
            FITS header containing all the infos of where the image is located (position, size, etc...)
        max_depth : int
            The moc resolution.
        mask : `numpy.ndarray`, optional
            A boolean array of the same size of the image where pixels having the value 1 are part of
            the final MOC and pixels having the value 0 are not.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            The MOC instance.
        """
        # Get the image dimensions
        height = header['NAXIS2']
        width = header['NAXIS1']

        # Get the image WCS to retrieve back its world coordinates
        wcs = wcs.WCS(header)

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
        world_pix = wcs.wcs_pix2world(pix, 1)

        # Call the from_lonlat method to create the MOC from the world coordinates.
        return MOCSkyRegion.from_lonlat(
            lon=world_pix[:, 0] * u.deg,
            lat=world_pix[:, 1] * u.deg,
            max_depth=max_depth,
            meta=meta,
            visual=visual)

    @classmethod
    def from_fits_images(cls, path_l, max_depth):
        """
        Loads a MOC from a set of FITS file images.

        Parameters
        ----------
        path_l : [str]
            A list of path where the fits image are located.
        max_depth : int
            The MOC resolution.

        Returns
        -------
        moc : `regions.MOCSkyRegion`
            The union of all the MOCs created from the paths found in ``path_l``.
        """
        moc = MOCSkyRegion()
        for path in path_l:
            header = fits.getheader(path)
            current_moc = MOCSkyRegion.from_image(header=header, max_depth=max_depth)
            moc = moc.union(current_moc)

        return moc

    @classmethod
    def from_skycoord(cls, skycoord, max_depth, meta=None, visual=None):
        """
        Creates a `~regions.MOCSkyRegion` from a `~astropy.coordinates.SkyCoord`.

        Parameters
        ----------
        skycoord : `astropy.coordinates.SkyCoord`
            The sky coordinates that will belong to the MOC.
        max_depth : int
            The depth of the smallest HEALPix cells contained in the MOC. Must be <= 29.
        meta : `~regions.RegionMeta` object, optional
            A dictionary which stores the meta attributes of this region.
        visual : `~regions.RegionVisual` object, optional
            A dictionary which stores the visual meta attributes of this region.

        Returns
        -------
        moc : `~regions.MOCSkyRegion`
            The MOC instance.
        """
        nside = level_to_nside(max_depth)
        hp = HEALPix(nside=(1 << max_depth), order='nested')
        ipix = hp.lonlat_to_healpix(skycoord.icrs.ra, skycoord.icrs.dec)

        shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - max_depth)
        itvs = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        itv = IntervalSet(itvs)
        return cls(itv, meta, visual)

    @classmethod
    def from_lonlat(cls, lon, lat, max_depth):
        """
        Creates a MOC from astropy lon, lat `astropy.units.Quantity`.

        Parameters
        ----------
        lon : `astropy.units.Quantity`
            The longitudes of the sky coordinates belonging to the MOC.
        lat : `astropy.units.Quantity`
            The latitudes of the sky coordinates belonging to the MOC.
        max_depth : int
            The depth of the smallest HEALPix cells contained in the MOC.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC
        """
        hp = HEALPix(nside=(1 << max_depth), order='nested')
        ipix = hp.lonlat_to_healpix(lon, lat)

        shift = 2 * (MOCSkyRegion.HPY_MAX_DEPTH - max_depth)
        itvs = np.vstack((ipix << shift, (ipix + 1) << shift)).T

        itv = IntervalSet(itvs)
        return cls(itv)

    @classmethod
    def from_polygon_skycoord(cls, skycoord, inside=None, max_depth=10):
        """
        Creates a MOC from a polygon.

        The polygon is given as an `astropy.coordinates.SkyCoord` that contains the
        vertices of the polygon. Concave and convex polygons are accepted but
        self-intersecting ones are currently not properly handled.

        Parameters
        ----------
        skycoord : `astropy.coordinates.SkyCoord`
            The sky coordinates defining the vertices of a polygon. It can describe a convex or
            concave polygon but not a self-intersecting one.
        inside : `astropy.coordinates.SkyCoord`, optional
            A point that will be inside the MOC is needed as it is not possible to determine the inside area of a polygon
            on the unit sphere (there is no infinite area that can be considered as the outside because on the sphere,
            a closed polygon delimits two finite areas).
            Possible improvement: take the inside area as the one covering the smallest region on the sphere.

            If inside=None (default behavior), the mean of all the vertices is taken as lying inside the polygon. That approach may not work for
            concave polygons.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC
        """
        return MOCSkyRegion.from_polygon(lon=skycoord.icrs.ra, lat=skycoord.icrs.dec,
                                inside=inside, max_depth=max_depth)

    @classmethod
    def from_polygon(cls, lon, lat, inside=None, max_depth=10):
        """
        Creates a MOC from a polygon

        The polygon is given as lon and lat `astropy.units.Quantity` that define the
        vertices of the polygon. Concave and convex polygons are accepted but
        self-intersecting ones are currently not properly handled.

        Parameters
        ----------
        lon : `astropy.units.Quantity`
            The longitudes defining the polygon. Can describe convex and
            concave polygons but not self-intersecting ones.
        lat : `astropy.units.Quantity`
            The latitudes defining the polygon. Can describe convex and concave
            polygons but not self-intersecting ones.
        inside : `astropy.coordinates.SkyCoord`, optional
            A point that will be inside the MOC is needed as it is not possible to determine the inside area of a polygon
            on the unit sphere (there is no infinite area that can be considered as the outside because on the sphere,
            a closed polygon delimits two finite areas).
            Possible improvement: take the inside area as the one covering the smallest region on the sphere.

            If inside=None (default behavior), the mean of all the vertices is taken as lying inside the polygon. That approach may not work for
            concave polygons.
        max_depth : int, optional
            The resolution of the MOC. Set to 10 by default.

        Returns
        -------
        result : `regions.MOCSkyRegion`
            The resulting MOC
        """
        from .polygon import PolygonComputer

        polygon_computer = PolygonComputer(lon, lat, inside, max_depth)
        # Create the moc from the python dictionary

        moc = MOCSkyRegion.from_json(polygon_computer.ipix)
        # We degrade it to the user-requested order
        if polygon_computer.degrade_to_max_depth:
            moc = moc.degrade_to_depth(max_depth)

        return moc

    @property
    def sky_fraction(self):
        """
        Sky fraction covered by the MOC
        """
        pix = self._best_res_pixels()
        num_pix = pix.size
        return num_pix / float(3 << (2*(self.max_depth + 1)))
