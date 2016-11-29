# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..point import PointPixelRegion, PointSkyRegion
from .utils import ASTROPY_LT_13


class TestPointPixelRegion:
    def test_str(self):
        center = PixCoord(3, 4)
        reg = PointPixelRegion(center)

        assert str(reg) == 'PointPixelRegion\ncenter: PixCoord(x=3, y=4)'


class TestPointSkyRegion:
    def test_str(self):
        center = SkyCoord(3, 4, unit='deg')
        reg = PointSkyRegion(center)

        if ASTROPY_LT_13:
            assert str(reg) == 'PointSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n    (3.0, 4.0)>'
        else:
            assert str(reg) == 'PointSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n    ( 3.,  4.)>'
