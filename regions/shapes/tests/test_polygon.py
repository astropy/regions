# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy import units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..polygon import PolygonPixelRegion, PolygonSkyRegion


def test_polygon_pixel():
    vertices = PixCoord([3, 4, 3], [3, 4, 4])
    reg = PolygonPixelRegion(vertices)

    assert str(reg) == 'PolygonPixelRegion\nvertices: PixCoord(x=[3 4 3], y=[3 4 4])'


def test_polygon_sky():
    vertices = SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg)
    reg = PolygonSkyRegion(vertices)

    expected = ('PolygonSkyRegion\nvertices: <SkyCoord (ICRS): (ra, dec) in deg\n'
                '    [(3.0, 3.0), (4.0, 4.0), (3.0, 4.0)]>')
    assert str(reg) == expected

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
