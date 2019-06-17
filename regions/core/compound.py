# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import operator as op

import numpy as np

from . import PixelRegion, SkyRegion, BoundingBox, RegionMask
from ..core.attributes import (CompoundRegionPix, CompoundRegionSky,
                               RegionMeta, RegionVisual)

__all__ = ['CompoundPixelRegion', 'CompoundSkyRegion']


class CompoundPixelRegion(PixelRegion):
    """
    Represents the logical combination of two regions in pixel coordinates.

    Parameters
    ----------
    region1 : `~regions.PixelRegion` object
        The inner Pixel region.
    region2 : `~regions.PixelRegion` object
        The outer Pixel region.
    operator : `function`
        A callable binary operator.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    region1 = CompoundRegionPix('region1')
    region2 = CompoundRegionPix('region2')

    def __init__(self, region1, region2, operator, meta=None, visual=None):

        if not callable(operator):
            raise TypeError("operator must be callable")

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
        self._repr_params = ('region1', 'region2', 'operator')

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
        bbox = BoundingBox(
            ixmin=min(mask1.bbox.ixmin, mask2.bbox.ixmin),
            ixmax=max(mask1.bbox.ixmax, mask2.bbox.ixmax),
            iymin=min(mask1.bbox.iymin, mask2.bbox.iymin),
            iymax=max(mask1.bbox.iymax, mask2.bbox.iymax)
        )

        # Pad mask1.data and mask2.data to get the same shape
        padded_data = list()
        for mask in (mask1, mask2):
            pleft = abs(mask.bbox.ixmin - bbox.ixmin)
            pright = abs(bbox.ixmax - mask.bbox.ixmax)
            ptop = abs(bbox.iymax - mask.bbox.iymax)
            pbottom = abs(mask.bbox.iymin - bbox.iymin)
            padded_data.append(np.pad(mask.data,
                                      ((ptop, pbottom), (pleft, pright)),
                                      'constant'))

        data = self.operator(*np.array(padded_data, dtype=np.int))
        return RegionMask(data=data, bbox=bbox)

    def to_sky(self, wcs):
        skyreg1 = self.region1.to_sky(wcs=wcs)
        skyreg2 = self.region2.to_sky(wcs=wcs)
        return CompoundSkyRegion(region1=skyreg1,
                                 operator=self.operator,
                                 region2=skyreg2, meta=self.meta, visual=self.visual)

    @staticmethod
    def _make_annulus_path(patch_inner, patch_outer):
        """
        Defines a matplotlib annulus path from two patches.

        This preserves the cubic Bezier curves (CURVE4) of the aperture
        paths.

        # This is borrowed from photutils aperture.
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
        Matplotlib patch object for annulus region (`matplotlib.patches.PathPatch`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.patches.PathPatch` object accepts

        Returns
        -------
        patch : `~matplotlib.patches.PathPatch`
            Matplotlib patch object
        """

        if self.region1.center == self.region2.center and self.operator == op.xor:
            import matplotlib.patches as mpatches

            patch_inner = self.region1.as_artist(origin=origin)
            patch_outer = self.region2.as_artist(origin=origin)
            path = self._make_annulus_path(patch_inner, patch_outer)
            patch = mpatches.PathPatch(path, **kwargs)
            return patch
        else:
            raise NotImplementedError

    def bounding_box(self, **kwargs):
        raise NotImplementedError

    @property
    def area(self):
        raise NotImplementedError

    def rotate(self, center, angle):
        """Make a rotated region.

        Rotates counter-clockwise for positive ``angle``.

        Parameters
        ----------
        center : `PixCoord`
            Rotation center point
        angle : `~astropy.coordinates.Angle`
            Rotation angle

        Returns
        -------
        region : `CompoundPixelRegion`
            Rotated region (an independent copy)
        """
        region1 = self.region1.rotate(center, angle)
        region2 = self.region2.rotate(center, angle)
        return self.copy(region1=region1, region2=region2)


class CompoundSkyRegion(SkyRegion):
    """
    Represents the logical combination of two regions in sky coordinates.

    Parameters
    ----------
    region1 : `~regions.SkyRegion` object
        The inner sky region.
    region2 : `~regions.SkyRegion` object
        The outer sky region.
    operator : `function`
        A callable binary operator.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    region1 = CompoundRegionSky('region1')
    region2 = CompoundRegionSky('region2')

    def __init__(self, region1, region2, operator, meta=None, visual=None):
        if not callable(operator):
            raise TypeError("operator must be callable")

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

        self._repr_params = ('region1', 'region2', 'operator')

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
        return CompoundPixelRegion(region1=pixreg1,
                                   operator=self.operator,
                                   region2=pixreg2, meta=self.meta, visual=self.visual)

    def as_artist(self, ax, **kwargs):
        raise NotImplementedError
