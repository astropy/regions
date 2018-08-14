# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..point import PointPixelRegion, PointSkyRegion
from .utils import ASTROPY_LT_13
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


class TestPointPixelRegion(BaseTestPixelRegion):

    reg = PointPixelRegion(PixCoord(3, 4))
    sample_box = [-2, 8, -1, 9]
    inside = []
    outside = [(3.1, 4.2), (5, 4)]
    expected_area = 0
    expected_repr = '<PointPixelRegion(PixCoord(x=3, y=4))>'
    expected_str = 'Region: PointPixelRegion\ncenter: PixCoord(x=3, y=4)'

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)


class TestPointSkyRegion(BaseTestSkyRegion):

    reg = PointSkyRegion(SkyCoord(3, 4, unit='deg'))

    if ASTROPY_LT_13:
        expected_repr = ('<PointSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                    '    (3.0, 4.0)>)>')
        expected_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                   '(ra, dec) in deg\n    (3.0, 4.0)>')
    else:
        expected_repr = ('<PointSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                    '    ( 3.,  4.)>)>')
        expected_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                   '(ra, dec) in deg\n    ( 3.,  4.)>')

    def test_contains(self):
        position = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        # points do not contain things
        assert reg.contains(position) == np.array([False,False], dtype='bool')
