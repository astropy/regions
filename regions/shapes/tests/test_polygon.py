# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_allclose, assert_equal
from astropy.tests.helper import pytest, assert_quantity_allclose
from astropy import units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from .utils import ASTROPY_LT_13


class TestPolygonPixelRegion:
    def setup(self):
        vertices = PixCoord([3, 4, 3], [3, 4, 4])
        self.poly = PolygonPixelRegion(vertices)

    def test_str(self):
        expected = 'PolygonPixelRegion\nvertices: PixCoord(x=[3 4 3], y=[3 4 4])'
        assert str(self.poly) == expected

    def _test_basic(self):
        """
        TODO: implement these tests!

        Just make sure that things don't crash, but no test of numerical accuracy
        """

        poly1 = PolygonPixelRegion(([10, 20, 20, 10, 10], [20, 20, 10, 10, 20]))

        poly2 = PolygonPixelRegion(([15, 25, 25, 15, 15], [23, 23, 13, 13, 23]))

        assert_allclose(poly1.area, 100)
        assert_allclose(poly2.area, 100)

        assert PixCoord(13, 12) in poly1
        assert not PixCoord(13, 12) in poly2

        coords = PixCoord([13, 8], [15, 15])
        assert_equal(poly1.contains(coords), [False, True])

        # poly3 = poly1.union(poly2)
        # poly4 = poly1.intersection(poly2)
        # assert_allclose(poly3.area, 165)
        # assert_allclose(poly4.area, 35)


class TestPolygonSkyRegion:
    def setup(self):
        vertices = SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg)
        self.poly = PolygonSkyRegion(vertices)

    def test_str(self):
        if ASTROPY_LT_13:
            expected = ('PolygonSkyRegion\nvertices: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    [(3.0, 3.0), (4.0, 4.0), (3.0, 4.0)]>')
        else:
            expected = ('PolygonSkyRegion\nvertices: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    [( 3.,  3.), ( 4.,  4.), ( 3.,  4.)]>')
        assert str(self.poly) == expected

    def _test_basic(self):
        """
        TODO: implement these tests!

        Just make sure that things don't crash, but no test of numerical accuracy
        """

        coords1 = SkyCoord([10, 20, 20, 10, 10] * u.deg, [20, 20, 10, 10, 20] * u.deg)
        poly1 = PolygonSkyRegion(coords1)

        coords1 = SkyCoord([15, 25, 25, 15, 15] * u.deg, [23, 23, 13, 13, 23] * u.deg)
        poly2 = PolygonSkyRegion(coords1)

        assert_quantity_allclose(poly1.area, 42 * u.sr)
        assert_quantity_allclose(poly2.area, 42 * u.sr)

        # TODO: test contains

        poly3 = poly1.union(poly2)
        poly4 = poly1.intersection(poly2)

        assert_quantity_allclose(poly3.area, 42 * u.sr)
        assert_quantity_allclose(poly4.area, 42 * u.sr)
