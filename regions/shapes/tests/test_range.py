# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import pytest
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose, assert_equal

from regions.core import RegionMeta, RegionVisual
from regions.core.compound import CompoundSphericalSkyRegion
from regions.shapes.circle import CircleSphericalSkyRegion
from regions.shapes.polygon import (PolygonPixelRegion, PolygonSkyRegion,
                                    PolygonSphericalSkyRegion)
from regions.shapes.range import RangeSphericalSkyRegion
from regions.shapes.tests.test_common import BaseTestSphericalSkyRegion


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestRangeSphericalSkyRegion(BaseTestSphericalSkyRegion):
    inside = [(5 * u.deg, 0 * u.deg)]
    outside = [(75 * u.deg, 10 * u.deg)]
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                  latitude_range=[-4, 4] * u.deg,
                                  frame='icrs',
                                  meta=meta, visual=visual)

    expected_repr = ('<RangeSphericalSkyRegion(\n'
                     'frame=icrs,\n'
                     'longitude_range=[ 0. 10.] deg,\n'
                     'latitude_range=[-4.  4.] deg\n'
                     ')>')
    expected_str = ('Region: RangeSphericalSkyRegion\n'
                    'frame: icrs\n'
                    'longitude_range: [ 0. 10.] deg\n'
                    'latitude_range: [-4.  4.] deg')

    def test_no_range_set(self):
        with pytest.raises(ValueError) as excinfo:
            _ = RangeSphericalSkyRegion(frame='icrs')
        estr = 'A range for at least one of longitude and latitude must be set'
        assert estr in str(excinfo.value)

    def test_lon_only(self):
        reg = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                      frame='icrs')
        assert reg._bound_nverts == 2
        assert reg.latitude_bounds is None

        coos = self.inside + self.outside
        lon, lat = zip(*coos)
        skycoord = SkyCoord(list(lon), list(lat))
        actual = reg.contains(skycoord)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        assert reg.contains(skycoord[0])

        reg_transf = reg.transform_to('galactic')
        skycoord_transf = skycoord.transform_to('galactic')
        actual = reg_transf.contains(skycoord_transf)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        # Manual "no latitude constraints"
        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[-90, 90] * u.deg,
                                       frame='icrs')
        assert reg2._bound_nverts == 2
        actual = reg2.contains(skycoord)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        assert reg2.contains(skycoord[0])

    def test_lat_only(self):
        reg = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                      frame='icrs')
        assert reg._bound_nverts == 0
        assert reg.longitude_bounds is None

        coos = self.inside + self.outside
        lon, lat = zip(*coos)
        skycoord = SkyCoord(list(lon), list(lat))
        actual = reg.contains(skycoord)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        assert reg.contains(skycoord[0])

        reg_transf = reg.transform_to('galactic')
        skycoord_transf = skycoord.transform_to('galactic')
        actual = reg_transf.contains(skycoord_transf)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

    def test_lon_over_meridian(self):
        reg = RangeSphericalSkyRegion(longitude_range=[358, 2] * u.deg,
                                      latitude_range=[-4, 4] * u.deg,
                                      frame='icrs')
        assert_allclose(reg.centroid.ra.deg, 0)
        assert_allclose(reg.centroid.dec.deg, 0)

    def test_lat_over_poles(self):
        reg = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                      latitude_range=[80, -80] * u.deg,
                                      frame='icrs')
        assert isinstance(reg.latitude_bounds, CompoundSphericalSkyRegion)
        assert reg._bound_nverts == 6

        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[90, -80] * u.deg,
                                       frame='icrs')
        assert isinstance(reg2.latitude_bounds, CircleSphericalSkyRegion)
        assert reg2.latitude_bounds.center.dec.deg == -90

        reg3 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[80, -90] * u.deg,
                                       frame='icrs')
        assert isinstance(reg3.latitude_bounds, CircleSphericalSkyRegion)
        assert reg3.latitude_bounds.center.dec.deg == 90

    def test_centroid(self):
        assert_allclose(self.reg.centroid.ra.deg, 5)
        assert_allclose(self.reg.centroid.dec.deg, 0)

        # Long and narrow range:
        reg = RangeSphericalSkyRegion(longitude_range=[0, 350] * u.deg,
                                      latitude_range=[45, 50] * u.deg,
                                      frame='icrs')
        assert reg.centroid == reg.centroid_avg

        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       frame='icrs')
        assert_allclose(reg2.centroid.ra.deg, 5)
        assert_allclose(reg2.centroid.dec.deg, 0)

        reg3 = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                       frame='icrs')
        assert_allclose(reg3.centroid.ra.deg, 0)
        assert_allclose(reg3.centroid.dec.deg, 90)

    def test_centroid_avg(self):
        assert_allclose(self.reg.centroid_avg.ra.deg, 5)
        assert_allclose(self.reg.centroid_avg.dec.deg, 0)

        # Very wide range:
        reg = RangeSphericalSkyRegion(longitude_range=[200, 100] * u.deg,
                                      latitude_range=[-50, 80] * u.deg,
                                      frame='icrs')
        assert reg.centroid_avg == SkyCoord(330 * u.deg, 15 * u.deg)

        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       frame='icrs')
        assert_allclose(reg2.centroid_avg.ra.deg, 5)
        assert_allclose(reg2.centroid_avg.dec.deg, 0)

        reg3 = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                       frame='icrs')
        assert_allclose(reg3.centroid_avg.ra.deg, 0)
        assert_allclose(reg3.centroid_avg.dec.deg, 90)

    def test_vertices(self):
        assert_allclose(self.reg.vertices.ra.deg, [0, 10, 10, 0])
        assert_allclose(self.reg.vertices.dec.deg, [-4, -4, 4, 4])

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.centroid.ra.deg, 5)
        assert_allclose(reg.centroid.dec.deg, 0)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        polypix = self.reg.to_pixel(wcs)
        assert isinstance(polypix, PolygonPixelRegion)
        assert len(polypix.vertices) == 4

        assert_allclose(polypix.vertices.x,
                        [-4536.82156523, -5685.54622545,
                         -5771.86292244, -4861.62586712])
        assert_allclose(polypix.vertices.y,
                        [-3091.16199865, -3236.05135267,
                         -2837.97224466, -2724.76777962])

        polysky = self.reg.to_sky(wcs)
        assert isinstance(polysky, PolygonSkyRegion)
        assert_allclose(self.reg.vertices.ra.deg, [0, 10, 10, 0])
        assert_allclose(self.reg.vertices.dec.deg, [-4, -4, 4, 4])
        assert len(polysky.vertices) == 4

        polypix2 = self.reg.to_pixel(wcs,
                                     include_boundary_distortions=True,
                                     discretize_kwargs={'n_points': 10})
        assert isinstance(polypix2, PolygonPixelRegion)
        assert len(polypix2.vertices) == 40

        polysky2 = self.reg.to_sky(wcs,
                                   include_boundary_distortions=True,
                                   discretize_kwargs={'n_points': 10})
        assert isinstance(polysky2, PolygonSkyRegion)
        assert len(polysky2.vertices) == 40

    def test_transformation_over_poles(self, wcs):
        # Transformation for over-poles case:
        reg = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                      latitude_range=[80, -80] * u.deg,
                                      frame='icrs')

        with pytest.raises(NotImplementedError):
            _ = reg.to_pixel(wcs)

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
        reg2 = self.reg.transform_to('galactic')
        reg2_cent = reg2.centroid
        transf_reg_cent = self.reg.centroid.transform_to('galactic')

        assert isinstance(reg2, RangeSphericalSkyRegion)
        assert not reg2._is_original_frame
        assert transf_reg_cent.frame.name == reg2_cent.frame.name
        assert_quantity_allclose(reg2_cent.l, transf_reg_cent.l)
        assert_quantity_allclose(reg2_cent.b, transf_reg_cent.b)
        assert reg2.frame.name == 'galactic'

        assert reg2 != self.reg

        reg3 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[80, 90] * u.deg,
                                       frame='icrs').transform_to('galactic')
        assert reg3 != reg2

        # Formally won't round trip, because the original lon/lat range info
        # is not stored:
        reg = reg2.transform_to('icrs')
        assert reg != self.reg

    def test_copy_frame_transformation(self):
        reg2 = self.reg.transform_to('galactic')
        reg3 = reg2.copy()

        assert reg2 == reg3

    def test_repr_str_frame_transformation(self):
        expected_repr_transf = ('<RangeSphericalSkyRegion(\n'
                                'frame=galactic,\n'
                                'longitude_bounds=<LuneSphericalSkyRegion('
                                'center_gc1=<SkyCoord (Galactic): (l, b) in deg\n'
                                '    (206.98913108, -11.42449097)>, '
                                'center_gc2=<SkyCoord (Galactic): (l, b) in deg\n'
                                '    (31.62719158, 2.54468138)>)>,\n'
                                'latitude_bounds=<CircleAnnulusSphericalSkyRegion('
                                'center=<SkyCoord (Galactic): (l, b) in deg\n'
                                '    (122.93192526, 27.12825241)>, '
                                'inner_radius=86.0 deg, outer_radius=94.0 deg)>\n'
                                ')>')
        expected_str_transf = ('Region: RangeSphericalSkyRegion\n'
                               'frame: galactic\n'
                               'longitude_bounds: <LuneSphericalSkyRegion('
                               'center_gc1=<SkyCoord (Galactic): (l, b) in deg\n'
                               '    (206.98913108, -11.42449097)>, '
                               'center_gc2=<SkyCoord (Galactic): (l, b) in deg\n'
                               '    (31.62719158, 2.54468138)>)>\n'
                               'latitude_bounds: <CircleAnnulusSphericalSkyRegion('
                               'center=<SkyCoord (Galactic): (l, b) in deg\n'
                               '    (122.93192526, 27.12825241)>, '
                               'inner_radius=86.0 deg, outer_radius=94.0 deg)>')

        reg2 = self.reg.transform_to('galactic')

        assert repr(reg2) == expected_repr_transf
        assert str(reg2) == expected_str_transf

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.longitude_range = [0, 4] * u.deg
        assert reg != self.reg

        reg2 = CircleSphericalSkyRegion(SkyCoord(5 * u.deg, 0 * u.deg),
                                        6.3999492595405645 * u.deg)
        assert self.reg != reg2

    def test_bounding_circle(self):
        skycoord = SkyCoord(5 * u.deg, 0 * u.deg)
        reg = CircleSphericalSkyRegion(skycoord,
                                       6.3999492595405645 * u.deg)

        bc = self.reg.bounding_circle
        assert_quantity_allclose(bc.radius, reg.radius)
        assert_quantity_allclose(bc.center.ra, reg.center.ra)
        assert_quantity_allclose(bc.center.dec, reg.center.dec)

        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       frame='icrs')
        assert reg2.bounding_circle == reg2.longitude_bounds.bounding_circle

        reg3 = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                       frame='icrs')
        assert reg3.bounding_circle == reg3.latitude_bounds.bounding_circle

    def test_bounding_circle_over_poles(self):
        # Wrap over poles:
        reg = RangeSphericalSkyRegion(longitude_range=[0, 350] * u.deg,
                                      latitude_range=[80, -80] * u.deg,
                                      frame='icrs')

        with pytest.raises(NotImplementedError):
            _ = reg.bounding_circle

    def test_bounding_lonlat(self):
        bounding_lonlat = self.reg.bounding_lonlat

        assert_quantity_allclose(bounding_lonlat[0],
                                 Longitude([0. * u.deg,
                                            10. * u.deg]))

        assert_quantity_allclose(bounding_lonlat[1],
                                 Latitude([-4 * u.deg,
                                           4 * u.deg]))

        # Transformed:
        reg_tr = self.reg.transform_to('galactic')
        bll_transf = reg_tr.bounding_lonlat
        assert_quantity_allclose(bll_transf[0],
                                 Longitude([92.7264313, 117.42725845] * u.deg))
        assert_quantity_allclose(bll_transf[1],
                                 Latitude([-66.71102705, -56.48535559] * u.deg))

        # Touching poles: should return longitude bounds
        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[80, 90] * u.deg,
                                       frame='icrs')
        bounding_lonlat2 = reg2.bounding_lonlat
        assert bounding_lonlat2[0] is not None

        reg3 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       latitude_range=[-90, -80] * u.deg,
                                       frame='icrs')
        bounding_lonlat3 = reg3.bounding_lonlat
        assert bounding_lonlat3[0] is not None

        # Only lon or lat ranges:
        reg4 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       frame='icrs')
        bll = reg4.bounding_lonlat
        bound_bll = reg4.longitude_bounds.bounding_lonlat
        assert_quantity_allclose(bll[0], bound_bll[0])
        assert_quantity_allclose(bll[1], bound_bll[1])

        reg5 = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                       frame='icrs')
        bll = reg5.bounding_lonlat
        bound_bll = reg5.latitude_bounds.bounding_lonlat
        assert bll[0] is None
        assert_quantity_allclose(bll[1], bound_bll[1])

    def test_bounding_lonlat_over_poles(self):
        # Wrap over poles:
        reg = RangeSphericalSkyRegion(longitude_range=[0, 350] * u.deg,
                                      latitude_range=[80, -80] * u.deg,
                                      frame='icrs')

        with pytest.raises(NotImplementedError):
            _ = reg.bounding_lonlat

    def test_discretize_boundary(self):
        polyrange = self.reg.discretize_boundary(n_points=100)
        assert isinstance(polyrange, PolygonSphericalSkyRegion)
        assert len(polyrange.vertices) == 400

        reg2 = RangeSphericalSkyRegion(longitude_range=[0, 10] * u.deg,
                                       frame='icrs')
        polyrange2 = reg2.discretize_boundary(n_points=100)
        assert polyrange2 == reg2.longitude_bounds.discretize_boundary(n_points=100)
        assert len(polyrange2.vertices) == 200

        reg3 = RangeSphericalSkyRegion(latitude_range=[-4, 4] * u.deg,
                                       frame='icrs')
        polyrange3 = reg3.discretize_boundary(n_points=100)
        assert polyrange3 == reg3.latitude_bounds.discretize_boundary(n_points=100)
        assert len(polyrange3.region1.vertices) == 100
        assert len(polyrange3.region2.vertices) == 100

    def test_discretize_over_poles(self):
        # Wrap over poles:
        reg = RangeSphericalSkyRegion(longitude_range=[0, 350] * u.deg,
                                      latitude_range=[80, -80] * u.deg,
                                      frame='icrs')

        with pytest.raises(NotImplementedError):
            _ = reg.discretize_boundary()
