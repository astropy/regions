# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.circle import (CirclePixelRegion, CircleSkyRegion,
                                   CircleSphericalSkyRegion)
from regions.shapes.polygon import (PolygonPixelRegion, PolygonSkyRegion,
                                    PolygonSphericalSkyRegion)
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion,
                                              BaseTestSphericalSkyRegion)
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

    def test_to_spherical_sky(self, wcs):
        sphskycircle = self.reg.to_spherical_sky(wcs,
                                                 include_boundary_distortions=False)
        assert isinstance(sphskycircle, CircleSphericalSkyRegion)

        try:
            sphskycircle = self.reg.to_spherical_sky(wcs,
                                                     include_boundary_distortions=True)
            assert isinstance(sphskycircle, PolygonSphericalSkyRegion)
        except NotImplementedError:
            pytest.xfail()

    def test_to_spherical_sky_no_wcs(self):
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_spherical_sky(include_boundary_distortions=True)
        estr = "'wcs' must be set if 'include_boundary_distortions'=True"
        assert estr in str(excinfo.value)

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
        assert_allclose(pixcircle.radius, 0.027777777777828305)

        skycircle2 = pixcircle.to_sky(wcs)

        assert_quantity_allclose(skycircle.center.data.lon,
                                 skycircle2.center.data.lon)
        assert_quantity_allclose(skycircle.center.data.lat,
                                 skycircle2.center.data.lat)
        assert_quantity_allclose(skycircle2.radius, skycircle.radius)

        sphskycircle = self.reg.to_spherical_sky(wcs,
                                                 include_boundary_distortions=False)
        assert isinstance(sphskycircle, CircleSphericalSkyRegion)

        try:
            sphskycircle = self.reg.to_spherical_sky(wcs,
                                                     include_boundary_distortions=True)
            assert isinstance(sphskycircle, PolygonSphericalSkyRegion)
        except NotImplementedError:
            pytest.xfail()

    def test_to_spherical_sky_no_wcs(self):
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_spherical_sky(include_boundary_distortions=True)
        estr = "'wcs' must be set if 'include_boundary_distortions'=True"
        assert estr in str(excinfo.value)

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


class TestCircleSphericalSkyRegion(BaseTestSphericalSkyRegion):
    inside = [(3 * u.deg, 4 * u.deg)]
    outside = [(3 * u.deg, 0 * u.deg)]
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CircleSphericalSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 2 * u.arcsec,
                                   meta=meta, visual=visual)

    expected_repr = ('<CircleSphericalSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    (3., 4.)>, radius=2.0 arcsec)>')
    expected_str = ('Region: CircleSphericalSkyRegion\ncenter: <SkyCoord (ICRS): '
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
        sphskycircle = CircleSphericalSkyRegion(skycoord, 2 * u.arcsec)

        pixcircle = sphskycircle.to_pixel(wcs)

        assert_allclose(pixcircle.center.x, -50.5)
        assert_allclose(pixcircle.center.y, 299.5)
        assert_allclose(pixcircle.radius, 0.027777777777828305)

        skycircle = sphskycircle.to_sky(wcs)
        assert isinstance(skycircle, CircleSkyRegion)

        sphskycircle2 = pixcircle.to_spherical_sky(wcs)

        assert_quantity_allclose(sphskycircle.center.data.lon,
                                 sphskycircle2.center.data.lon)
        assert_quantity_allclose(sphskycircle.center.data.lat,
                                 sphskycircle2.center.data.lat)
        assert_quantity_allclose(sphskycircle2.radius, sphskycircle.radius)

        polysky = sphskycircle.to_sky(wcs, include_boundary_distortions=True)
        assert isinstance(polysky, PolygonSkyRegion)

        polypix = sphskycircle.to_pixel(wcs, include_boundary_distortions=True)
        assert isinstance(polypix, PolygonPixelRegion)

    def test_transformation_no_wcs(self):
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_sky(include_boundary_distortions=True)
        estr = "'wcs' must be set if 'include_boundary_distortions'=True"
        assert estr in str(excinfo.value)

        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_pixel(include_boundary_distortions=True)
        estr = "'wcs' must be set if 'include_boundary_distortions'=True"
        assert estr in str(excinfo.value)

    def test_frame_transformation(self):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        reg = CircleSphericalSkyRegion(skycoord, 2 * u.arcsec)

        reg2 = reg.transform_to('icrs')
        assert reg2.center == skycoord.transform_to('icrs')
        assert_allclose(reg2.radius.to_value('arcsec'), 2)
        assert isinstance(reg2, CircleSphericalSkyRegion)
        assert reg2.frame.name == 'icrs'

    def test_dimension_center(self):
        center = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        radius = 2 * u.arcsec
        with pytest.raises(ValueError) as excinfo:
            CircleSphericalSkyRegion(center, radius)
        estr = "'center' must be a scalar SkyCoord"
        assert estr in str(excinfo.value)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.radius = 3 * u.arcsec
        assert reg != self.reg

    def test_zero_size(self):
        with pytest.raises(ValueError):
            CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 0. * u.arcsec)

    def test_bounding_circle(self):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        reg = CircleSphericalSkyRegion(skycoord, 2 * u.arcsec)

        bounding_circle = reg.bounding_circle
        assert bounding_circle == reg

    def test_bounding_lonlat(self):
        skycoord = SkyCoord(3 * u.deg, 0 * u.deg, frame='galactic')
        reg = CircleSphericalSkyRegion(skycoord, 2 * u.arcsec)
        bounding_lonlat = reg.bounding_lonlat

        assert_quantity_allclose(bounding_lonlat[0],
                                 Longitude([3. * u.deg - 2 * u.arcsec,
                                            3. * u.deg + 2 * u.arcsec]))

        assert_quantity_allclose(bounding_lonlat[1],
                                 Latitude([-2 * u.arcsec,
                                           2 * u.arcsec]))
