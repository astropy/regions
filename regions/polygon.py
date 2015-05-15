from astropy import units as u
from astropy.coordinates import SkyCoord

from spherical_geometry.polygon import SphericalPolygon

from .core import SkyRegion


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
