# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_allclose
from astropy.tests.helper import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB


class TestEllipsePixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = EllipsePixelRegion(
            center=center,
            major=4,
            minor=3,
            angle=5 * u.deg,
        )

    def test_str(self):
        expected = 'EllipsePixelRegion\ncenter: PixCoord(x=3, y=4)\nmajor: 4\nminor: 3\nangle: 5.0 deg'
        assert str(self.reg) == expected

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.width, 8)
        assert_allclose(patch.height, 6)
        assert_allclose(patch.angle, 5)


class TestEllipseSkyRegion:
    def setup(self):
        center = SkyCoord(3, 4, unit='deg')
        self.reg = EllipseSkyRegion(
            center=center,
            major=4 * u.deg,
            minor=3 * u.deg,
            angle=5 * u.deg,
        )

    def test_str(self):

        if ASTROPY_LT_13:
            expected = ('EllipseSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    (3.0, 4.0)>\nmajor: 4.0 deg\nminor: 3.0 deg\nangle: 5.0 deg')
        else:
            expected = ('EllipseSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    ( 3.,  4.)>\nmajor: 4.0 deg\nminor: 3.0 deg\nangle: 5.0 deg')

        assert str(self.reg) == expected
