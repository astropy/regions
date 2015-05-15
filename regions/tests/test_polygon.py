from astropy import units as u
from astropy.coordinates import SkyCoord

from ..polygon import SkyPolygonRegion


def test_basic():
    """
    Just make sure that things don't crash, but no test of numerical accuracy
    """

    coords1 = SkyCoord([10, 20, 20, 10, 10] * u.deg, [20, 20, 10, 10, 20] * u.deg)
    poly1 = SkyPolygonRegion(coords1)

    coords1 = SkyCoord([15, 25, 25, 15, 15] * u.deg, [23, 23, 13, 13, 23] * u.deg)
    poly2 = SkyPolygonRegion(coords1)

    poly3 = poly1.union(poly2)
    poly4 = poly1.intersection(poly2)

    poly1.area
    poly2.area
    poly3.area
    poly4.area
