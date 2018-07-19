# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from . import PixelRegion, SkyRegion, BoundingBox, Mask
from ..core.attributes import CompoundRegionPix, CompoundRegionSky


__all__ = ['CompoundPixelRegion', 'CompoundSkyRegion']


class CompoundPixelRegion(PixelRegion):
    """
    Represents the logical combination of two regions in pixel coordinates.
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
            self.meta = meta
        if visual is None:
            self.visual = region1.visual
        else:
            self.visual = visual
        self._operator = operator
        self._repr_params = ('region1', 'region2', 'operator')

    @property
    def operator(self):
        return self._operator

    def contains(self, pixcoord):
        in_reg = self.operator(self.region1.contains(pixcoord), self.region2.contains(pixcoord))
        if self.meta.get('inverted', False):
            return not in_reg
        else:
            return in_reg

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
        return Mask(data=data, bbox=bbox)

    def to_sky(self, wcs):
        skyreg1 = self.region1.to_sky(wcs=wcs)
        skyreg2 = self.region2.to_sky(wcs=wcs)
        return CompoundSkyRegion(region1=skyreg1,
                                 operator=self.operator,
                                 region2=skyreg2, meta=self.meta, visual=self.visual)

    def as_patch(self, **kwargs):
        raise NotImplementedError

    def to_shapely(self, **kwargs):
        raise NotImplementedError

    def bounding_box(self, **kwargs):
        raise NotImplementedError

    @property
    def area(self):
        raise NotImplementedError


class CompoundSkyRegion(SkyRegion):
    """
    Represents the logical combination of two regions in sky coordinates.
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
            self.meta = meta
        if visual is None:
            self.visual = region1.visual
        else:
            self.visual = visual
        self._operator = operator

        self._repr_params = ('region1', 'region2', 'operator')

    @property
    def operator(self):
        return self._operator

    def contains(self, skycoord, wcs):
        in_reg = self.operator(self.region1.contains(skycoord, wcs),
                             self.region2.contains(skycoord, wcs))
        if self.meta.get('inverted', False):
            return not in_reg
        else:
            return in_reg

    def to_pixel(self, wcs):
        pixreg1 = self.region1.to_pixel(wcs=wcs)
        pixreg2 = self.region2.to_pixel(wcs=wcs)
        return CompoundPixelRegion(region1=pixreg1,
                                   operator=self.operator,
                                   region2=pixreg2, meta=self.meta, visual=self.visual)

    def as_patch(self, ax, **kwargs):
        raise NotImplementedError
