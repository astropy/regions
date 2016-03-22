from astropy import units as u
from astropy.coordinates import SkyCoord

from spherical_geometry.polygon import SphericalPolygon

from .core import SkyRegion


import math

import numpy as np
from astropy import units as u

from .core import PixelRegion


class PolygonRegion(object):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    pixcoord : tuple of `numpy.ndarray`
        The positions of the polygon vertices, as a tuple of two
        `numpy.ndarray` instances. In future this could also be a `PixCoord`
        instance.
    """

    def __init__(self, pixcoord):
        self.x, self.y = pixcoord

    @property
    def area(self):
        return self.to_shapely().area

    def __contains__(self, pixcoord):
        from shapely.geometry import Point
        x, y = pixcoord
        shape = self.to_shapely()
        if np.isscalar(x):
            return shape.contains(Point(x, y))
        else:
            return np.array([shape.contains(Point(x[i], y[i])) for i in range(len(x))])

    def to_shapely(self):
        from shapely.geometry import Polygon
        return Polygon(zip(self.x, self.y))

    @classmethod
    def from_shapely(cls, polygon):
        # TODO: currently we assume polygons have no holes
        return cls(zip(*polygon.exterior.coords))

    def union(self, other):
        shape1 = self.to_shapely()
        shape2 = other.to_shapely()
        return PolygonRegion.from_shapely(shape1.union(shape2))

    def intersection(self, other):
        shape1 = self.to_shapely()
        shape2 = other.to_shapely()
        return PolygonRegion.from_shapely(shape1.intersection(shape2))

    def symmetric_difference(self, other):
        shape1 = self.to_shapely()
        shape2 = other.to_shapely()
        return PolygonRegion.from_shapely(shape1.symmetric_difference(shape2))


def _skycoord_to_spherical_polygon(skycoord):
    return SphericalPolygon.from_radec(skycoord.spherical.lon.degree,
                                       skycoord.spherical.lat.degree,
                                       degrees=True)

def _spherical_polygon_to_skycoord(polygon, frame):
    lon, lat = tuple(polygon.to_radec())[0]
    return SkyCoord(lon, lat, unit=u.degree, frame=frame)


class SkyPolygonRegion(SkyRegion):
    """
    A spherical polygon on the sky.

    Parameters
    ----------
    skycoord : `~astropy.coordinates.SkyCoord`
        The coordinates of the vertices of the polygon
    """

    def __init__(self, skycoord):
        self.frame = skycoord.frame
        self.polygon = _skycoord_to_spherical_polygon(skycoord)

    @property
    def area(self):
        return self.polygon.area() * u.sr

    def union(self, other):
        result = self.polygon.union(other.polygon)
        return SkyPolygonRegion(_spherical_polygon_to_skycoord(result, self.frame))

    def intersection(self, other):
        result = self.polygon.intersection(other.polygon)
        return SkyPolygonRegion(_spherical_polygon_to_skycoord(result, self.frame))
