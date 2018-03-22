# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from . import PixCoord, PixelRegion, SkyRegion, BoundingBox, Mask
import numpy as np

__all__ = ['CompoundPixelRegion', 'CompoundSkyRegion']


class CompoundPixelRegion(PixelRegion):
    """
    Represents the logical combination of two regions in pixel coordinates.
    """

    def __init__(self, region1, region2, operator):
        if not isinstance(region1, PixelRegion):
            raise TypeError("region1 must be a PixelRegion")
        if not isinstance(region2, PixelRegion):
            raise TypeError("region2 must be a PixelRegion")
        if not callable(operator):
            raise TypeError("operator must be callable")

        self.region1 = region1
        self.region2 = region2
        self.operator = operator
        self._repr_params = [('component 1', self.region1),
                             ('component 2', self.region2),
                             ('operator', self.operator),
                            ]

    def contains(self, pixcoord):
        raise NotImplementedError

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

        bbox_borders = np.array([bbox.ixmin, bbox.ixmax, bbox.iymin, bbox.iymax]) 
        if (bbox_borders < 0).any():
            raise NotImplementedError("Bounding box must be within array for "
                                      "compound regions, see ")

        # Pad mask1.data and mask2.data to get the same shape
        padded_data = list()
        for mask in (mask1, mask2):
            pleft = mask.bbox.ixmin - bbox.ixmin
            pright = bbox.ixmax - mask.bbox.ixmax
            ptop = bbox.iymax - mask.bbox.iymax
            pbottom = mask.bbox.iymin - bbox.iymin
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
                                 region2=skyreg2)
        raise NotImplementedError

    def as_patch(self, **kwargs):
        raise NotImplementedError

    def to_shapely(self, **kwargs):
        raise NotImplementedError

    def bounding_box(self, **kwargs):
        raise NotImplementedError


class CompoundSkyRegion(SkyRegion):
    """
    Represents the logical combination of two regions in sky coordinates.
    """

    def __init__(self, region1, region2, operator):
        if not isinstance(region1, SkyRegion):
            raise TypeError("region1 must be a SkyRegion")
        if not isinstance(region2, SkyRegion):
            raise TypeError("region2 must be a SkyRegion")
        if not callable(operator):
            raise TypeError("operator must be callable")

        self.region1 = region1
        self.region2 = region2
        self.operator = operator

        self._repr_params = [('component 1', self.region1),
                             ('component 2', self.region2),
                             ('operator', self.operator),
                            ]

    def contains(self, skycoord, wcs):
        return self.operator(self.region1.contains(skycoord, wcs),
                             self.region2.contains(skycoord, wcs))

    def to_pixel(self, wcs):
        pixreg1 = self.region1.to_pixel(wcs=wcs)
        pixreg2 = self.region2.to_pixel(wcs=wcs)
        return CompoundPixelRegion(region1=pixreg1,
                                   operator=self.operator,
                                   region2=pixreg2)

    def as_patch(self, ax, **kwargs):
        raise NotImplementedError

    @property
    def area(self):
        raise NotImplementedError
