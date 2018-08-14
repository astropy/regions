# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..text import TextPixelRegion, TextSkyRegion
from .utils import ASTROPY_LT_13
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


class TestTextPixelRegion(BaseTestPixelRegion):

    reg = TextPixelRegion(PixCoord(3, 4), "Sample Text")
    sample_box = [-2, 8, -1, 9]
    inside = []
    outside = [(3.1, 4.2), (5, 4)]
    expected_area = 0
    expected_repr = '<TextPixelRegion(PixCoord(x=3, y=4), text=Sample Text)>'
    expected_str = 'Region: TextPixelRegion\ncenter: PixCoord(x=3, y=4)\ntext: Sample Text'

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)


class TestTextSkyRegion(BaseTestSkyRegion):

    reg = TextSkyRegion(SkyCoord(3, 4, unit='deg'), "Sample Text")

    if ASTROPY_LT_13:
        expected_repr = ('<TextSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                            '    (3.0, 4.0)>, text=Sample Text)>')
        expected_str = ('Region: TextSkyRegion\ncenter: <SkyCoord (ICRS): '
                          '(ra, dec) in deg\n    (3.0, 4.0)>\ntext: Sample Text')
    else:
        expected_repr = ('<TextSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                             '    ( 3.,  4.)>, text=Sample Text)>')
        expected_str = ('Region: TextSkyRegion\ncenter: <SkyCoord (ICRS): '
                        '(ra, dec) in deg\n    ( 3.,  4.)>\ntext: Sample Text')

    def test_contains(self, wcs):
        position = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        # Nothing is in a text region
        assert all(self.reg.contains(position, wcs) == np.array([False, False], dtype='bool'))
