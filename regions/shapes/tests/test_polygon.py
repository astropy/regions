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
        # We will be using this polygon for basic tests:
        #
        # (3,0) *
        #       |\
        #       | \
        #       |  \
        # (0,0) *---* (2, 0)
        #
        vertices = PixCoord([0, 2, 0], [0, 0, 3])
        self.reg = PolygonPixelRegion(vertices)
        self.pixcoord_inside = PixCoord(1, 1)
        self.pixcoord_outside = PixCoord(1, 2)

    def test_repr_str(self):
        reg_repr = ('<PolygonPixelRegion(vertices=PixCoord(x=[0 2 0], '
                    'y=[0 0 3]))>')
        assert repr(self.reg) == reg_repr

        reg_str = ('Region: PolygonPixelRegion\nvertices: PixCoord(x=[0 2 0],'
                   ' y=[0 0 3])')
        assert str(self.reg) == reg_str

    def test_contains_scalar(self):
        assert self.reg.contains(self.pixcoord_inside)
        assert self.pixcoord_inside in self.reg

        assert not self.reg.contains(self.pixcoord_outside)
        assert self.pixcoord_outside not in self.reg

    def test_contains_array_1d(self):
        pixcoord = PixCoord([1, 1], [1, 2])
        actual = self.reg.contains(pixcoord)
        expected = [True, False]
        assert_equal(actual, expected)

        with pytest.raises(ValueError) as exc:
            pixcoord in self.reg
        assert 'coord must be scalar' in str(exc)

    # TODO: fix implementation to work with 2D arrays
    @pytest.mark.xfail
    def test_contains_array_2d(self):
        pixcoord = PixCoord(
            [[1, 1, 1], [1, 1, 1]],
            [[1, 1, 1], [2, 2, 2]],
        )
        actual = self.reg.contains(pixcoord)
        expected = [[True, True, True], [False, False, False]]
        assert_equal(actual, expected)

    def _test_basic(self):
        poly1 = PolygonPixelRegion(([10, 20, 20, 10, 10], [20, 20, 10, 10, 20]))

        poly2 = PolygonPixelRegion(([15, 25, 25, 15, 15], [23, 23, 13, 13, 23]))

        # assert_allclose(poly1.area, 100)
        # assert_allclose(poly2.area, 100)

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

    def test_repr_str(self):
        if ASTROPY_LT_13:
            reg_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                        'dec) in deg\n    [(3.0, 3.0), (4.0, 4.0), (3.0, '
                        '4.0)]>)>')
            reg_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
                       ' (ra, dec) in deg\n    [(3.0, 3.0), (4.0, 4.0), (3.0,'
                       ' 4.0)]>')
        else:
            reg_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                        'dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,  '
                        '4.)]>)>')
            reg_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
                       ' (ra, dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,'
                       '  4.)]>')

        assert repr(self.poly) == reg_repr
        assert str(self.poly) == reg_str

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
