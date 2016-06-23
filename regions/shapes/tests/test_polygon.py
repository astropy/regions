import math
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from ...core import PixCoord
from ..polygon import PolygonPixelRegion, PolygonSkyRegion


def test_init_pixel():
    vertices = PixCoord([3, 4, 3], [3, 4, 4])
    c = PolygonPixelRegion(vertices)


def test_init_sky():
    vertices = SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg)
    c = PolygonSkyRegion(vertices)

# def test_basic():
#     """
#     Just make sure that things don't crash, but no test of numerical accuracy
#     """
#
#     poly1 = PolygonRegion(([10, 20, 20, 10, 10], [20, 20, 10, 10, 20]))
#
#     poly2 = PolygonRegion(([15, 25, 25, 15, 15], [23, 23, 13, 13, 23]))
#
#     poly3 = poly1.union(poly2)
#     poly4 = poly1.intersection(poly2)
#
#     np.testing.assert_allclose(poly1.area, 100)
#     np.testing.assert_allclose(poly2.area, 100)
#     np.testing.assert_allclose(poly3.area, 165)
#     np.testing.assert_allclose(poly4.area, 35)
#
#     assert (13, 12) in poly1
#     assert not (13, 12) in poly2
#
#     # TODO: it seems __contains__ has to return a scalar value?
#     # assert np.all(((np.array([13, 8]), np.array([15, 15])) in poly1) == np.array([False, True]))
#
#
# def test_sky_basic():
#     """
#     Just make sure that things don't crash, but no test of numerical accuracy
#     """
#
#     coords1 = SkyCoord([10, 20, 20, 10, 10] * u.deg, [20, 20, 10, 10, 20] * u.deg)
#     poly1 = SkyPolygonRegion(coords1)
#
#     coords1 = SkyCoord([15, 25, 25, 15, 15] * u.deg, [23, 23, 13, 13, 23] * u.deg)
#     poly2 = SkyPolygonRegion(coords1)
#
#     poly3 = poly1.union(poly2)
#     poly4 = poly1.intersection(poly2)
#
#     poly1.area
#     poly2.area
#     poly3.area
#     poly4.area
