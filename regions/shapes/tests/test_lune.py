# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import pytest
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions.core import RegionMeta, RegionVisual
from regions.shapes.circle import CircleSphericalSkyRegion
from regions.shapes.lune import LuneSphericalSkyRegion
from regions.shapes.polygon import PolygonSphericalSkyRegion
from regions.shapes.tests.test_common import BaseTestSphericalSkyRegion


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestLuneSphericalSkyRegion(BaseTestSphericalSkyRegion):
    inside = [(90 * u.deg, 0 * u.deg)]
    outside = [(75 * u.deg, 0 * u.deg)]
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = LuneSphericalSkyRegion(SkyCoord(3 * u.deg, 0 * u.deg),
                                 SkyCoord(178 * u.deg, 0 * u.deg),
                                 meta=meta, visual=visual)

    expected_repr = ('<LuneSphericalSkyRegion(center_gc1=<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    (3., 0.)>, center_gc2=<SkyCoord (ICRS): (ra, dec) in deg\n'
                     '    (178., 0.)>)>')
    expected_str = ('Region: LuneSphericalSkyRegion\ncenter_gc1: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 0.)>\ncenter_gc2: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (178., 0.)>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center_gc1.ra.deg, 3)
        assert_allclose(reg.center_gc2.ra.deg, 178)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        # Test error for no boundary distortions:
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_pixel(wcs)
        estr = 'Invalid parameter: `include_boundary_distortions=False`!'
        assert estr in str(excinfo.value)

        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_sky(wcs)
        estr = 'Invalid parameter: `include_boundary_distortions=False`!'
        assert estr in str(excinfo.value)

        # Test for not implemented with include_boundary_distortions=True
        with pytest.raises(NotImplementedError):
            _ = self.reg.to_sky(wcs,
                                include_boundary_distortions=True,
                                discretize_kwargs={'n_points': 2})

        with pytest.raises(NotImplementedError):
            _ = self.reg.to_pixel(wcs,
                                  include_boundary_distortions=True,
                                  discretize_kwargs={'n_points': 2})

    def test_frame_transformation(self):
        reg2 = self.reg.transform_to('galactic')
        reg2_cent = reg2.centroid
        transf_reg_cent = self.reg.centroid.transform_to('galactic')
        assert isinstance(reg2, LuneSphericalSkyRegion)
        assert transf_reg_cent.frame.name == reg2_cent.frame.name
        assert_quantity_allclose(reg2_cent.l, transf_reg_cent.l)
        assert_quantity_allclose(reg2_cent.b, transf_reg_cent.b)
        assert_allclose(reg2.center_gc1.l.deg,
                        self.reg.center_gc1.transform_to('galactic').l.deg)
        assert reg2.frame.name == 'galactic'

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.center_gc1 = SkyCoord(4 * u.deg, 0 * u.deg)
        assert reg != self.reg

    def test_bounding_circle(self):
        skycoord = SkyCoord(90.5 * u.deg, 0 * u.deg)
        reg = CircleSphericalSkyRegion(skycoord, 90 * u.deg)

        bc = self.reg.bounding_circle
        assert bc == reg

    def test_bounding_lonlat(self):
        bounding_lonlat = self.reg.bounding_lonlat

        assert_quantity_allclose(bounding_lonlat[0],
                                 Longitude([88. * u.deg,
                                            93. * u.deg]))

        assert_quantity_allclose(bounding_lonlat[1],
                                 Latitude([-90 * u.deg,
                                           90 * u.deg]))

        reg2 = LuneSphericalSkyRegion(SkyCoord(3 * u.deg, 20 * u.deg),
                                      SkyCoord(178 * u.deg, 20 * u.deg),
                                      meta=self.meta, visual=self.visual)
        bounding_lonlat2 = reg2.bounding_lonlat
        assert bounding_lonlat2[0] is None

        reg3 = LuneSphericalSkyRegion(SkyCoord(3 * u.deg, -20 * u.deg),
                                      SkyCoord(178 * u.deg, -20 * u.deg),
                                      meta=self.meta, visual=self.visual)
        bounding_lonlat3 = reg3.bounding_lonlat
        assert bounding_lonlat3[0] is None

    def test_discretize_boundary(self):
        polylune = self.reg.discretize_boundary(n_points=100)
        assert isinstance(polylune, PolygonSphericalSkyRegion)
        assert len(polylune.vertices) == 200
