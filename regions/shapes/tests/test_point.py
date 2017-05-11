# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..point import PointPixelRegion, PointSkyRegion
from .utils import ASTROPY_LT_13
from ...tests.helpers import make_simple_wcs
from numpy.testing import assert_allclose
from astropy import units as u


class TestPointPixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = PointPixelRegion(center)

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)

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
