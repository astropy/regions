# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_allclose
from astropy.tests.helper import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB


class TestRectanglePixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = RectanglePixelRegion(
            center=center,
            width=4,
            height=3,
            angle=5 * u.deg,
        )

    def test_str(self):
        expected = 'RectanglePixelRegion\ncenter: PixCoord(x=3, y=4)\nwidth: 4\nheight: 3\nangle: 5.0 deg'
        assert str(self.reg) == expected

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.xy, (3, 4))
        assert_allclose(patch.get_width(), 4)
        assert_allclose(patch.get_height(), 3)
        assert_allclose(patch._angle, 5)


class TestRectangleSkyRegion:
    def setup(self):
        center = SkyCoord(3, 4, unit='deg')
        self.reg = RectangleSkyRegion(
            center=center,
            width=4 * u.deg,
            height=3 * u.deg,
            angle=5 * u.deg,
        )

    def test_str(self):
        if ASTROPY_LT_13:
            expected = ('RectangleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    (3.0, 4.0)>\nwidth: 4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')
        else:
            expected = ('RectangleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    ( 3.,  4.)>\nwidth: 4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')
        assert str(self.reg) == expected
