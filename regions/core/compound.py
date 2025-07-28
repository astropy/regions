# Licensed under a 3-clause BSD style license - see LICENSE.rst
import operator as op

import numpy as np

from regions.core.attributes import RegionType
from regions.core.core import PixelRegion, SkyRegion, SphericalSkyRegion
from regions.core.mask import RegionMask

__all__ = ['CompoundPixelRegion', 'CompoundSkyRegion', 'CompoundSphericalSkyRegion']


class CompoundPixelRegion(PixelRegion):
    """
    A class that represents the logical combination of two regions in
    pixel coordinates.

    Parameters
    ----------
    region1 : `~regions.PixelRegion`
        The inner Pixel region.
    region2 : `~regions.PixelRegion`
        The outer Pixel region.
    operator : callable
        A callable binary operator.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('region1', 'region2', 'operator')
    _mpl_artist = 'Patch'
    region1 = RegionType('region1', PixelRegion)
    region2 = RegionType('region2', PixelRegion)

    def __init__(self, region1, region2, operator, meta=None, visual=None):
        if not callable(operator):
            raise TypeError('operator must be callable')

        self.region1 = region1
        self.region2 = region2
        if meta is None:
            self.meta = region1.meta
        else:
            self.meta = meta
        if visual is None:
            self.visual = region1.visual
        else:
            self.visual = visual
        self._operator = operator

    @property
    def operator(self):
        return self._operator

    def contains(self, pixcoord):
        in_reg = self.operator(self.region1.contains(pixcoord),
                               self.region2.contains(pixcoord))
        if self.meta.get('include', True):
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_mask(self, mode='center', subpixels=1):
        if mode != 'center':
            raise NotImplementedError

        mask1 = self.region1.to_mask(mode=mode, subpixels=subpixels)
        mask2 = self.region2.to_mask(mode=mode, subpixels=subpixels)

        # Common bounding box
        bbox = self.bounding_box

        # Pad mask1.data and mask2.data to get the same shape
        padded_data = list()
        for mask in (mask1, mask2):
            pleft = abs(mask.bbox.ixmin - bbox.ixmin)
            pright = abs(bbox.ixmax - mask.bbox.ixmax)
            ptop = abs(bbox.iymax - mask.bbox.iymax)
            pbottom = abs(mask.bbox.iymin - bbox.iymin)
            padded_data.append(np.pad(mask.data,
                                      ((pbottom, ptop), (pleft, pright)),
                                      'constant'))

        data = self.operator(*np.array(padded_data, dtype=int))
        return RegionMask(data=data, bbox=bbox)

    def to_sky(self, wcs):
        skyreg1 = self.region1.to_sky(wcs=wcs)
        skyreg2 = self.region2.to_sky(wcs=wcs)
        return CompoundSkyRegion(region1=skyreg1, operator=self.operator,
                                 region2=skyreg2, meta=self.meta.copy(),
                                 visual=self.visual.copy())

    @staticmethod
    def _make_annulus_path(patch_inner, patch_outer):
        """
        Define a matplotlib annulus path from two patches.

        This preserves the cubic Bezier curves (CURVE4) of the aperture
        paths.

        Taken from ``photutils.aperture.core``.
        """
        import matplotlib.path as mpath

        path_inner = patch_inner.get_path()
        transform_inner = patch_inner.get_transform()
        path_inner = transform_inner.transform_path(path_inner)

        path_outer = patch_outer.get_path()
        transform_outer = patch_outer.get_transform()
        path_outer = transform_outer.transform_path(path_outer)

        verts_inner = path_inner.vertices[:-1][::-1]
        verts_inner = np.concatenate((verts_inner, [verts_inner[-1]]))

        verts = np.vstack((path_outer.vertices, verts_inner))
        codes = np.hstack((path_outer.codes, path_inner.codes))

        return mpath.Path(verts, codes)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.PathPatch`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : `dict`
            Any keyword arguments accepted by
            `~matplotlib.patches.PathPatch`.

        Returns
        -------
        patch : `~matplotlib.patches.PathPatch`
            A matplotlib patch object.
        """
        if (self.region1.center == self.region2.center
                and self.operator is op.xor):

            import matplotlib.patches as mpatches

            # set mpl_kwargs before as_artist is called on region1 and
            # region2
            mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
            mpl_kwargs.update(kwargs)

            patch_inner = self.region1.as_artist(origin=origin)
            patch_outer = self.region2.as_artist(origin=origin)
            path = self._make_annulus_path(patch_inner, patch_outer)

            patch = mpatches.PathPatch(path, **mpl_kwargs)
            return patch
        else:
            raise ValueError('unable to convert region to matplotlib '
                             f'artist: {self}')

    @property
    def bounding_box(self):
        return self.region1.bounding_box | self.region2.bounding_box

    @property
    def area(self):
        raise NotImplementedError

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `CompoundPixelRegion`
            The rotated region (which is an independent copy).
        """
        region1 = self.region1.rotate(center, angle)
        region2 = self.region2.rotate(center, angle)
        return self.copy(region1=region1, region2=region2)


class CompoundSkyRegion(SkyRegion):
    """
    A class that represents the logical combination of two regions in
    sky coordinates.

    Parameters
    ----------
    region1 : `~regions.SkyRegion`
        The inner sky region.
    region2 : `~regions.SkyRegion`
        The outer sky region.
    operator : callable
        A callable binary operator.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('region1', 'region2', 'operator')
    region1 = RegionType('region1', SkyRegion)
    region2 = RegionType('region2', SkyRegion)

    def __init__(self, region1, region2, operator, meta=None, visual=None):
        if not callable(operator):
            raise TypeError('operator must be callable')

        self.region1 = region1
        self.region2 = region2
        if meta is None:
            self.meta = region1.meta
        else:
            self.meta = meta
        if visual is None:
            self.visual = region1.visual
        else:
            self.visual = visual
        self._operator = operator

    @property
    def operator(self):
        return self._operator

    def contains(self, skycoord, wcs):
        in_reg = self.operator(self.region1.contains(skycoord, wcs),
                               self.region2.contains(skycoord, wcs))
        if self.meta.get('include', True):
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_pixel(self, wcs):
        pixreg1 = self.region1.to_pixel(wcs=wcs)
        pixreg2 = self.region2.to_pixel(wcs=wcs)
        return CompoundPixelRegion(region1=pixreg1, operator=self.operator,
                                   region2=pixreg2, meta=self.meta.copy(),
                                   visual=self.visual.copy())

    def as_artist(self, ax, **kwargs):
        raise NotImplementedError


class CompoundSphericalSkyRegion(SphericalSkyRegion):
    """
    A class that represents the logical combination of two regions in
    spherical sky coordinates.

    Parameters
    ----------
    region1 : `~regions.SphericalSkyRegion`
       First spherical sky region.
    region2 : `~regions.SphericalSkyRegion`
        Second spherical sky region.
    operator : callable
        A callable binary operator.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('region1', 'region2', 'operator')
    region1 = RegionType('region1', SphericalSkyRegion)
    region2 = RegionType('region2', SphericalSkyRegion)

    def __init__(self, region1, region2, operator, meta=None, visual=None):
        if not callable(operator):
            raise TypeError('operator must be callable')

        self.region1 = region1
        self.region2 = region2
        if meta is None:
            self.meta = region1.meta
        else:
            self.meta = RegionMeta()
        if visual is None:
            self.visual = region1.visual
        else:
            self.visual = RegionVisual()
        self._operator = operator

    @property
    def operator(self):
        return self._operator

    @property
    def frame(self):
        return self.region1.frame

    @property
    def bounding_circle(self):
        from regions.shapes.circle import CircleSphericalSkyRegion

        if self.operator in [op.or_]:
            # Union:
            # Compute connection between constituent two bounding circles,
            # then total "span" on either side, then half and get center.

            circ1 = self.region1.bounding_circle
            circ2 = self.region2.bounding_circle

            pa = circ1.center.position_angle(circ2.center)

            e1 = circ1.center.directional_offset_by(180 + pa, circ1.radius)
            e2 = circ2.center.directional_offset_by(pa, circ2.radius)
            sep = e1.separation(e2)

            cnew = e1.directional_offset_by(pa, sep / 2)

            return CircleSphericalSkyRegion(cnew, sep / 2)

        # Disjoint set / XOR:
        # Not clean
        # Intersection:
        # NOT anywhere near as clean. Much better to define for each subclass
        # that has a compound region as their internal region logic.
        # TODO:
        # General solution requires finding intersections of boundaries,
        # and then making a bounding circle. Doable with composed
        # circle boundaries & polygon-specific logic from spherical-geometry
        # Maybe better to not make bounding_circle and bounding_lonlat
        # abstract classes of the base SphericalSkyRegion, and only
        # add as properties to individual explicitly defined classes?
        raise NotImplementedError

    def contains(self, coord):
        in_reg = self.operator(
            self.region1.contains(coord), self.region2.contains(coord)
        )
        if self.meta.get('include', True):
            return in_reg
        else:
            return np.logical_not(in_reg)

    def transform_to(self, frame, merge_attributes=True):
        # Transform the consitituent regions then recompose:
        return self.__class__(
            self.region1.transform_to(frame, merge_attributes=merge_attributes),
            self.region2.transform_to(frame, merge_attributes=merge_attributes),
            self.operator,
            meta=self.meta,
            visual=self.visual,
        )

    def discretize_boundary(self, n_points=100):
        return CompoundSphericalSkyRegion(
            self.region1.discretize_boundary(n_points=n_points),
            self.region2.discretize_boundary(n_points=n_points),
            operator=self.operator,
            meta=self.meta.copy(),
            visual=self.visual.copy()
        )

    def to_sky(
        self, wcs=None, include_boundary_distortions=False,
        discretize_kwargs=None,
    ):
        planarreg1 = self.region1.to_sky(
            wcs=wcs,
            include_boundary_distortions=include_boundary_distortions,
            discretize_kwargs=discretize_kwargs,
        )
        planarreg2 = self.region2.to_sky(
            wcs=wcs,
            include_boundary_distortions=include_boundary_distortions,
            discretize_kwargs=discretize_kwargs,
        )
        return CompoundSkyRegion(
            region1=planarreg1,
            region2=planarreg2,
            operator=self.operator,
            meta=self.meta.copy(),
            visual=self.visual.copy(),
        )

    def to_pixel(
        self, wcs=None, include_boundary_distortions=False,
        discretize_kwargs=None,
    ):
        pixreg1 = self.region1.to_pixel(
            wcs=wcs,
            include_boundary_distortions=include_boundary_distortions,
            discretize_kwargs=discretize_kwargs,
        )
        pixreg2 = self.region2.to_pixel(
            wcs=wcs,
            include_boundary_distortions=include_boundary_distortions,
            discretize_kwargs=discretize_kwargs,
        )
        return CompoundPixelRegion(
            region1=pixreg1,
            region2=pixreg2,
            operator=self.operator,
            meta=self.meta.copy(),
            visual=self.visual.copy(),
        )
