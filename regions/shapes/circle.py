import math

import numpy as np
from astropy import units as u

from .core import PixelRegion


class CircularRegion(object):
    """
    A circle or set of circles in pixel coordinates.

    Parameters
    ----------
    pixcoord : tuple
        The position or positions of the center of the circle, as a tuple of
        scalars or arrays. In future this could also be a `PixCoord` instance.
    radius : float
        The radius of the circles
    """

    def __init__(self, pixcoord, radius):
        self.x, self.y = pixcoord
        self.radius = radius

    @property
    def area(self):
        return math.pi * self.radius ** 2

    @property
    def isscalar(self):
        return np.isscalar(self.x)

    def __contains__(self, pixcoord):
        # TODO: this should not work if both pixcoord and the regions are
        # arrays
        x, y = pixcoord
        return np.hypot(x - self.x, y - self.y) < self.radius

    def to_shapely(self):
        if not self.isscalar:
            raise TypeError("Can only convert single regions to shapely objects")
        from shapely.geometry import Point
        return Point(self.x, self.y).buffer(self.radius)

