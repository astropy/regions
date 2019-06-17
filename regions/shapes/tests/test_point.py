# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..point import PointPixelRegion, PointSkyRegion
from .utils import HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


@pytest.fixture(scope='session')
def wcs():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestPointPixelRegion(BaseTestPixelRegion):

    reg = PointPixelRegion(PixCoord(3, 4))
    sample_box = [-2, 8, -1, 9]
    inside = []
    outside = [(3.1, 4.2), (5, 4)]
    expected_area = 0
    expected_repr = '<PointPixelRegion(PixCoord(x=3, y=4))>'
    expected_str = 'Region: PointPixelRegion\ncenter: PixCoord(x=3, y=4)'

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_artist(self):
        artist = self.reg.as_artist()

        assert artist.get_data() == ([3], [4])

        artist = self.reg.as_artist(origin=(1, 1))

        assert artist.get_data() == ([2], [3])

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))


class TestPointSkyRegion(BaseTestSkyRegion):

    reg = PointSkyRegion(SkyCoord(3, 4, unit='deg'))

    expected_repr = ('<PointSkyRegion(<SkyCoord (ICRS): (ra, dec) in deg\n'
                     '    ( 3.,  4.)>)>')
    expected_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    ( 3.,  4.)>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_contains(self, wcs):
        position = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        # points do not contain things
        assert all(self.reg.contains(position, wcs) == np.array([False,False], dtype='bool'))
