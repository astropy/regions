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

    def test_contains_array_2d(self):
        pixcoord = PixCoord(
            [[1, 1, 1], [1, 1, 1]],
            [[1, 1, 1], [2, 2, 2]],
        )
        actual = self.reg.contains(pixcoord)
        expected = [[True, True, True], [False, False, False]]
        assert_equal(actual, expected)


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
