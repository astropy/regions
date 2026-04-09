# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.circle import CirclePixelRegion, CircleSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


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

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
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

    def test_zero_size(self):
        with pytest.raises(ValueError):
            CirclePixelRegion(PixCoord(50, 50), radius=0)


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
        assert_allclose(reg.radius.to_value('arcsec'), 2)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        skycircle = CircleSkyRegion(skycoord, 2 * u.arcsec)

        pixcircle = skycircle.to_pixel(wcs)

        assert_allclose(pixcircle.center.x, -50.5)
        assert_allclose(pixcircle.center.y, 299.5)
        assert_allclose(pixcircle.radius, 0.0278117, rtol=1e-6)

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
        estr = "'center' must be a scalar SkyCoord"
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

    def test_zero_size(self):
        with pytest.raises(ValueError):
            CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 0. * u.arcsec)


class TestCirclePixelRegionToSkyEllipse:
    """
    Tests for CirclePixelRegion.to_sky with use_ellipse=True.
    """

    def setup_method(self):
        self.center = PixCoord(10, 10)
        self.radius = 5.0
        self.meta = RegionMeta({'text': 'test'})
        self.visual = RegionVisual({'color': 'blue'})
        self.reg = CirclePixelRegion(self.center, self.radius,
                                     meta=self.meta, visual=self.visual)
        self.wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg),
                                   0.1 * u.deg, 20)

    def test_to_sky_use_ellipse(self):
        result = self.reg.to_sky(self.wcs, use_ellipse=True)
        assert isinstance(result, EllipseSkyRegion)
        assert result.meta == self.meta
        assert result.visual == self.visual

    def test_to_sky_ellipse_roundtrip(self):
        sky_ellipse = self.reg.to_sky(self.wcs, use_ellipse=True)
        pix_ellipse = sky_ellipse.to_pixel(self.wcs)
        # For a simple WCS without distortion, the roundtrip should
        # recover the original diameter as width and height
        assert_allclose(pix_ellipse.width, 2 * self.radius, rtol=1e-4)
        assert_allclose(pix_ellipse.height, 2 * self.radius, rtol=1e-4)
        assert_allclose(pix_ellipse.center.x, self.center.x, rtol=1e-5)
        assert_allclose(pix_ellipse.center.y, self.center.y, rtol=1e-5)

    def test_to_sky_ellipse_center_matches_circle(self):
        sky_circle = self.reg.to_sky(self.wcs)
        sky_ellipse = self.reg.to_sky(self.wcs, use_ellipse=True)
        assert_quantity_allclose(sky_ellipse.center.ra,
                                 sky_circle.center.ra)
        assert_quantity_allclose(sky_ellipse.center.dec,
                                 sky_circle.center.dec)

    def test_to_sky_use_ellipse_meta_copies(self):
        result = self.reg.to_sky(self.wcs, use_ellipse=True)
        result.meta['text'] = 'new'
        result.visual['color'] = 'green'
        assert result.meta['text'] != self.reg.meta['text']
        assert result.visual['color'] != self.reg.visual['color']


class TestCircleSkyRegionToPixelEllipse:
    """
    Tests for CircleSkyRegion.to_pixel with use_ellipse=True.
    """

    def setup_method(self):
        self.center = SkyCoord(3 * u.deg, 4 * u.deg)
        self.radius = 20 * u.arcsec
        self.meta = RegionMeta({'text': 'test'})
        self.visual = RegionVisual({'color': 'blue'})
        self.reg = CircleSkyRegion(self.center, self.radius,
                                   meta=self.meta, visual=self.visual)
        self.wcs = make_simple_wcs(SkyCoord(3 * u.deg, 4 * u.deg),
                                   5 * u.arcsec, 20)

    def test_to_pixel_use_ellipse(self):
        result = self.reg.to_pixel(self.wcs, use_ellipse=True)
        assert isinstance(result, EllipsePixelRegion)
        assert result.meta == self.meta
        assert result.visual == self.visual

    def test_to_pixel_default_returns_circle(self):
        result = self.reg.to_pixel(self.wcs)
        assert isinstance(result, CirclePixelRegion)

    def test_to_pixel_ellipse_roundtrip(self):
        pix_ellipse = self.reg.to_pixel(self.wcs, use_ellipse=True)
        sky_ellipse = pix_ellipse.to_sky(self.wcs)
        # Roundtrip should recover the original diameter
        assert_quantity_allclose(sky_ellipse.width,
                                 2 * self.radius, rtol=1e-4)
        assert_quantity_allclose(sky_ellipse.height,
                                 2 * self.radius, rtol=1e-4)

    def test_to_pixel_ellipse_center_matches_circle(self):
        pix_circle = self.reg.to_pixel(self.wcs)
        pix_ellipse = self.reg.to_pixel(self.wcs, use_ellipse=True)
        assert_allclose(pix_ellipse.center.x, pix_circle.center.x)
        assert_allclose(pix_ellipse.center.y, pix_circle.center.y)

    def test_to_pixel_use_ellipse_meta_copies(self):
        result = self.reg.to_pixel(self.wcs, use_ellipse=True)
        result.meta['text'] = 'new'
        result.visual['color'] = 'green'
        assert result.meta['text'] != self.reg.meta['text']
        assert result.visual['color'] != self.reg.visual['color']
