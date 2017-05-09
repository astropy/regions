# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from .core import PixelRegion, SkyRegion

__all__ = ['CompoundPixelRegion', 'CompoundSkyRegion']


class CompoundPixelRegion(PixelRegion):
    """
    Represents the logical combination of two regions in pixel coordinates.
    """

    def __init__(self, region1, operator, region2):
        self.region1 = region1
        self.region2 = region2
        self.operator = operator
        self._repr_params = [('component 1', self.region1),
                             ('component 2', self.region2),
                             ('operator', self.operator),
                            ]

    def contains(self, pixcoord):
        raise NotImplementedError

    def to_mask(self, mode='center'):
        raise NotImplementedError

    def to_sky(self, wcs, mode='local', tolerance=None):
        raise NotImplementedError

    def as_patch(self, **kwargs):
        raise NotImplementedError

    def to_shapely(self, **kwargs):
        raise NotImplementedError

    def __repr__(self):
        return "({0} {1} {2})".format(self.region1, self.operator, self.region2)


class CompoundSkyRegion(SkyRegion):
    """
    Represents the logical combination of two regions in sky coordinates.
    """

    def __init__(self, region1, region2, operator):
        self.region1 = region1
        self.region2 = region2
        self.operator = operator
        self._repr_params = [('component 1', self.region1),
                             ('component 2', self.region2),
                             ('operator', self.operator),
                            ]

    def contains(self, skycoord):
        return self.operator(self.region1.contains(skycoord),
                             self.region2.contains(skycoord))

    def to_pixel(self, wcs, mode='local', tolerance=None):
        raise NotImplementedError

    def as_patch(self, ax, **kwargs):
        raise NotImplementedError

    def __repr__(self):
        return "({0}\n{1}\n{2})".format(self.region1, self.operator, self.region2)

    @property
    def area(self):
        raise NotImplementedError
