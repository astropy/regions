# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from ..core import PixCoord, PixelRegion, SkyRegion, BoundingBox

__all__ = ['PointPixelRegion', 'PointSkyRegion', 'TextSkyRegion', 'TextPixelRegion']


class PointPixelRegion(PixelRegion):
    """
    A point position in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the point
    """

    def __init__(self, center, meta=None, visual=None):
        self.center = PixCoord._validate(center, name='center', expected='scalar')
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = None

    @property
    def area(self):
        return 0

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            return False
        else:
            return np.zeros(pixcoord.x.shape, dtype=bool)

    def to_shapely(self):
        return self.center.to_shapely()

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return PointSkyRegion(center)

    @property
    def bounding_box(self):
        return BoundingBox.from_float(self.center.x, self.center.x,
                                      self.center.y, self.center.y)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_patch(self, **kwargs):
        # FIXME: need to make radius constant
        from matplotlib.patches import Circle
        return Circle((self.center.x, self.center.y), radius=2, **kwargs)


class PointSkyRegion(SkyRegion):
    """
    A pixel region in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the point
    """

    def __init__(self, center, meta=None, visual=None):
        if center.isscalar:
            self.center = center
        else:
            raise ValueError('the centre should be a 0D SkyCoord object')
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = None

    def contains(self, skycoord, wcs):
        return False

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return PointPixelRegion(center)
