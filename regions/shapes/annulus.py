# Licensed under a 3-clause BSD style license - see LICENSE.rst
import math
import operator

import numpy as np

from astropy import wcs
from astropy import coordinates
from astropy import units as u

from ..core import CompoundPixelRegion, CompoundSkyRegion
from ..shapes import CirclePixelRegion, CircleSkyRegion

__all__ = ['CircleAnnulusPixelRegion', 'CircleAnnulusSkyRegion']


class CircleAnnulusPixelRegion(CompoundPixelRegion):
    """
    A circular annulus in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the annulus.
    inner_radius : float
        The inner radius of the annulus
    outer_radius : float
        The outer radius of the annulus
    """

    def __init__(self, center, inner_radius, outer_radius, meta=None, visual=None):
        region1 = CirclePixelRegion(center, inner_radius)
        region2 = CirclePixelRegion(center, outer_radius)
        super(CircleAnnulusPixelRegion, self).__init__(
            region1=region1, operator=operator.xor, region2=region2)
        self._repr_params = [('inner radius', region1.radius),
                             ('outer radius', region2.radius),
                             ('center', region2.center)]


    @property
    def center(self):
        return self.region1.center

    @property
    def inner_radius(self):
        return self.region1.radius

    @property
    def outer_radius(self):
        return self.region2.radius

    def bounding_box(self):
        return self.region2.bounding_box()

    def to_shapely(self):
        r1 = self.region1.to_shapely()
        r2 = self.region2.to_shapely()
        return r2.symmetric_difference(r1)


class CircleAnnulusSkyRegion(CompoundSkyRegion):
    """
    A circular annulus in sky coordinates.

    Parameters
    ----------
    center : :class:`~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_radius : :class:`~astropy.units.Quantity`
        The inner radius of the annulus in angular units
    outer_radius : :class:`~astropy.units.Quantity`
        The outer radius of the annulus in angular units
    """

    def __init__(self, center, inner_radius, outer_radius, meta=None, visual=None):
        region1 = CircleSkyRegion(center, inner_radius)
        region2 = CircleSkyRegion(center, outer_radius)
        super(CircleAnnulusSkyRegion, self).__init__(
            region1=region1, operator=operator.xor, region2=region2)
        self._repr_params = [('inner radius', region1.radius),
                             ('outer radius', region2.radius),
                             ('center', region2.center)]

    @property
    def center(self):
        return self.region1.center

    @property
    def inner_radius(self):
        return self.region1.radius

    @property
    def outer_radius(self):
        return self.region2.radius
