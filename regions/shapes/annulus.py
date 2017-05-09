# Licensed under a 3-clause BSD style license - see LICENSE.rst
import math
import operator

import numpy as np

from astropy import wcs
from astropy import coordinates
from astropy import units as u

from ..core import CompoundPixelRegion, CompoundSkyRegion
from ..shapes import CirclePixelRegion, CircleSkyRegion


class AnnulusPixelRegion(CompoundPixelRegion):
    """
    An annulus in pixel coordinates.

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
        super(AnnulusPixelRegion, self).__init__(
            region1, region2, operator.xor)
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

    def bounding_box():
        return self.region2.bounding_box()


class AnnulusSkyRegion(CompoundSkyRegion):
    """
    An annulus in sky coordinates.

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
        super(AnnulusSkyRegion, self).__init__(
            region1, region2, operator.xor)
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
