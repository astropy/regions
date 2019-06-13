# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from .utils import HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


@pytest.fixture(scope='session')
def wcs():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestEllipsePixelRegion(BaseTestPixelRegion):
    reg = EllipsePixelRegion(center=PixCoord(3, 4), width=4, height=3, angle=5 * u.deg)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 4)]
    expected_area = 3 * np.pi
    expected_repr = '<EllipsePixelRegion(PixCoord(x=3, y=4), width=4, height=3, angle=5.0 deg)>'
    expected_str = ('Region: EllipsePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
                    'width: 4\nheight: 3\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.width == 4
        assert reg.height == 3
        assert_allclose(reg.angle.to_value("deg"), 5)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.width, self.reg.width)
        assert_allclose(reg_new.height, self.reg.height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.width, 4)
        assert_allclose(patch.height, 3)
        assert_allclose(patch.angle, 5)

    def test_region_bbox(self):
        a = 7
        b = 3
        reg = EllipsePixelRegion(PixCoord(50, 50), width=a, height=b,
                                 angle=0. * u.deg)
        assert reg.bounding_box.shape == (b, a)

        reg = EllipsePixelRegion(PixCoord(50.5, 50.5), width=a, height=b,
                                 angle=0. * u.deg)
        assert reg.bounding_box.shape == (b + 1, a + 1)

        reg = EllipsePixelRegion(PixCoord(50, 50), width=a, height=b,
                                 angle=90. * u.deg)
        assert reg.bounding_box.shape == (a, b)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value("deg"), 95)


class TestEllipseSkyRegion(BaseTestSkyRegion):
    reg = EllipseSkyRegion(
        center=SkyCoord(3, 4, unit='deg'),
        width=4 * u.deg,
        height=3 * u.deg,
        angle=5 * u.deg,
    )

    expected_repr = ('<EllipseSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    ( 3.,  4.)>, width=4.0 deg, height=3.0 deg,'
                     ' angle=5.0 deg)>')
    expected_str = ('Region: EllipseSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    ( 3.,  4.)>\nwidth: 4.0 deg\n'
                    'height: 3.0 deg\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.width.to_value("deg"), 4)
        assert_allclose(reg.height.to_value("deg"), 3)
        assert_allclose(reg.angle.to_value("deg"), 5)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_dimension_center(self):
        center = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        width = 2 * u.arcsec
        height = 3 * u.arcsec
        with pytest.raises(ValueError) as err:
            EllipseSkyRegion(center, width, height)
        assert 'The center must be a 0D SkyCoord object' in str(err)

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs) == np.array([False, True], dtype='bool'))
