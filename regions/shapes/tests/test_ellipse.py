# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from .utils import ASTROPY_LT_13


class TestEllipsePixelRegion:
    def test_str(self):
        center = PixCoord(3, 4)
        reg = EllipsePixelRegion(center, 3, 4, 5 * u.deg)

        assert str(reg) == 'EllipsePixelRegion\ncenter: PixCoord(x=3, y=4)\nminor: 3\nmajor: 4\nangle: 5.0 deg'


class TestEllipseSkyRegion:
    def test_str(self):
        center = SkyCoord(3, 4, unit='deg')
        reg = EllipseSkyRegion(center, 3 * u.deg, 4 * u.deg, 5 * u.deg)

        if ASTROPY_LT_13:
            expected = ('EllipseSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    (3.0, 4.0)>\nminor: 3.0 deg\nmajor: 4.0 deg\nangle: 5.0 deg')
        else:
            expected = ('EllipseSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    ( 3.,  4.)>\nminor: 3.0 deg\nmajor: 4.0 deg\nangle: 5.0 deg')

        assert str(reg) == expected
