import math

import numpy as np
from astropy import units as u

from ..core import PixelRegion, SkyRegion


class CirclePixelRegion(PixelRegion):
    """
    A circle in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the circle.
    radius : float
        The radius of the circle
    """

    def __init__(self, center, radius,):
        self.center = center
        self.radius = radius

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def __contains__(self, pixcoord):
        return np.hypot(x - self.center.x, y - self.center.y) < self.radius

    def to_shapely(self):
        return self.center.to_shapely().buffer(self.radius)

    def to_sky(self, wcs, mode='local', tolerance=None):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def to_mask(self, mode='center'):
        # TOOD: needs to be implemented
        raise NotImplementedError("")



class CircleSkyRegion(SkyRegion):
    """
    A circle in sky coordinates.

    Parameters
    ----------
    center : :class:`~astropy.coordinates.SkyCoord`
        The position of the center of the circle.
    radius : :class:`~astropy.units.Quantity`
        The radius of the circle in angular units
    """

    def __init__(self, center, radius):
        self.center = center
        self.radius = radius
    
    @property
    def area(self):
        return math.pi * self.radius ** 2

    def __contains__(self, skycoord):
        return self.center.separation(skycoord)

    def to_pixel(self, wcs, mode='local', tolerance=None):
        # TOOD: needs to be implemented
        raise NotImplementedError("")
