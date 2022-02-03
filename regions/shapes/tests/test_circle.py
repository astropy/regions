# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import numpy as np
from numpy.testing import assert_allclose
import pytest

from ...core import PixCoord, RegionMeta, RegionVisual
from ...tests.helpers import make_simple_wcs
from ..circle import CirclePixelRegion, CircleSkyRegion
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion
from ..utils import HAS_MATPLOTLIB  # noqa


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestCirclePixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CirclePixelRegion(PixCoord(3, 4), radius=2, meta=meta, visual=visual)
    sample_box = [0, 6, 1, 7]
    inside = [(3, 4)]
    outside = [(3, 0)]
    expected_area = 4 * np.pi
    expected_repr = '<CirclePixelRegion(center=PixCoord(x=3, y=4), radius=2)>'
    expected_str = ('Region: CirclePixelRegion\ncenter: PixCoord(x=3, y=4)'
                    '\nradius: 2')

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.radius == 2
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.radius, self.reg.radius)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.radius, 2)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.radius = 3
        assert reg != self.reg


class TestCircleSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 2 * u.arcsec,
                          meta=meta, visual=visual)

    expected_repr = ('<CircleSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    (3., 4.)>, radius=2.0 arcsec)>')
    expected_str = ('Region: CircleSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>\nradius: 2.0 '
                    'arcsec')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.radius.to_value("arcsec"), 2)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        skycircle = CircleSkyRegion(skycoord, 2 * u.arcsec)

        pixcircle = skycircle.to_pixel(wcs)

        assert_allclose(pixcircle.center.x, -50.5)
        assert_allclose(pixcircle.center.y, 299.5)
        assert_allclose(pixcircle.radius, 0.027777777777828305)

        skycircle2 = pixcircle.to_sky(wcs)

        assert_quantity_allclose(skycircle.center.data.lon,
                                 skycircle2.center.data.lon)
        assert_quantity_allclose(skycircle.center.data.lat,
                                 skycircle2.center.data.lat)
        assert_quantity_allclose(skycircle2.radius, skycircle.radius)

    def test_dimension_center(self):
        center = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        radius = 2 * u.arcsec
        with pytest.raises(ValueError) as excinfo:
            CircleSkyRegion(center, radius)
        estr = 'The center must be a scalar SkyCoord object'
        assert estr in str(excinfo.value)

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, True], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.radius = 3 * u.arcsec
        assert reg != self.reg
