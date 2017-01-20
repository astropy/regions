# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..point import PointPixelRegion, PointSkyRegion
from .utils import ASTROPY_LT_13


class TestPointPixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = PointPixelRegion(center)

    def test_repr_str(self):
        reg_repr = '<PointPixelRegion(PixCoord(x=3, y=4))>'
        assert repr(self.reg) == reg_repr

        reg_str = 'Region: PointPixelRegion\ncenter: PixCoord(x=3, y=4)'
        assert str(self.reg) == reg_str


class TestPointSkyRegion:
    def setup(self):
        center = SkyCoord(3, 4, unit='deg')
        self.reg = PointSkyRegion(center)

    def test_repr_str(self):
        if ASTROPY_LT_13:
            reg_repr = ('<PointSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    (3.0, 4.0)>)>')
            reg_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    (3.0, 4.0)>')
        else:
            reg_repr = ('<PointSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                        '    ( 3.,  4.)>)>')
            reg_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    ( 3.,  4.)>')

        assert repr(self.reg) == reg_repr
        assert str(self.reg) == reg_str
