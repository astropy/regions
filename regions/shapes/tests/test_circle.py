# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose, assert_equal

from regions._utils.optional_deps import HAS_MATPLOTLIB, HAS_SHAPELY
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.circle import (CirclePixelRegion, CircleSkyRegion,
                                   CircleSphericalSkyRegion)
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.polygon import PolygonPixelRegion, PolygonSkyRegion
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

        with pytest.raises(NotImplementedError):
            _ = self.reg.to_spherical_sky(wcs,
                                          include_boundary_distortions=True)

    def test_to_spherical_sky_no_wcs(self):
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_spherical_sky(include_boundary_distortions=True)
        estr = "'wcs' must be set if `include_boundary_distortions=True`"
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


def test_contains_covers_ignore_include_metadata():
    """
    Test that contains and covers ignore the 'include' metadata key, and
    depend strictly on the geometric shape.
    """
    meta = RegionMeta({'include': False})
    circle_exclude = CirclePixelRegion(PixCoord(10, 10), radius=5, meta=meta)
    circle_include = CirclePixelRegion(PixCoord(10, 10), radius=5)

    inside_pt = PixCoord(10, 10)
    boundary_pt = PixCoord(15, 10)
    outside_pt = PixCoord(20, 10)

    # Should be completely identical regardless of 'include' metadata
    assert_equal(circle_exclude.contains(inside_pt),
                 circle_include.contains(inside_pt))
    assert_equal(circle_exclude.contains(boundary_pt),
                 circle_include.contains(boundary_pt))
    assert_equal(circle_exclude.contains(outside_pt),
                 circle_include.contains(outside_pt))

    # Should be completely identical regardless of 'include' metadata
    assert_equal(circle_exclude.covers(inside_pt),
                 circle_include.covers(inside_pt))
    assert_equal(circle_exclude.covers(boundary_pt),
                 circle_include.covers(boundary_pt))
    assert_equal(circle_exclude.covers(outside_pt),
                 circle_include.covers(outside_pt))


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

        sphskycircle = self.reg.to_spherical_sky(wcs,
                                                 include_boundary_distortions=False)
        assert isinstance(sphskycircle, CircleSphericalSkyRegion)

        with pytest.raises(NotImplementedError):
            _ = self.reg.to_spherical_sky(wcs,
                                          include_boundary_distortions=True)

    def test_to_spherical_sky_no_wcs(self):
        with pytest.raises(ValueError) as excinfo:
            _ = self.reg.to_spherical_sky(include_boundary_distortions=True)
        estr = "'wcs' must be set if `include_boundary_distortions=True`"
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
        assert_allclose(pixcircle.radius, 0.0278117, rtol=1e-6)

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
        estr = "'wcs' must be set if `include_boundary_distortions=True`"
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


class TestCirclePixelRegionToSkyEllipse:
    """
    Tests for CirclePixelRegion.to_sky with as_ellipse=True.
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

    def test_to_sky_as_ellipse(self):
        result = self.reg.to_sky(self.wcs, as_ellipse=True)
        assert isinstance(result, EllipseSkyRegion)
        assert result.meta == self.meta
        assert result.visual == self.visual

    def test_to_sky_ellipse_roundtrip(self):
        sky_ellipse = self.reg.to_sky(self.wcs, as_ellipse=True)
        pix_ellipse = sky_ellipse.to_pixel(self.wcs)
        # For a simple WCS without distortion, the roundtrip should
        # recover the original diameter as width and height
        assert_allclose(pix_ellipse.width, 2 * self.radius, rtol=1e-4)
        assert_allclose(pix_ellipse.height, 2 * self.radius, rtol=1e-4)
        assert_allclose(pix_ellipse.center.x, self.center.x, rtol=1e-5)
        assert_allclose(pix_ellipse.center.y, self.center.y, rtol=1e-5)

    def test_to_sky_ellipse_center_matches_circle(self):
        sky_circle = self.reg.to_sky(self.wcs)
        sky_ellipse = self.reg.to_sky(self.wcs, as_ellipse=True)
        assert_quantity_allclose(sky_ellipse.center.ra,
                                 sky_circle.center.ra)
        assert_quantity_allclose(sky_ellipse.center.dec,
                                 sky_circle.center.dec)

    def test_to_sky_as_ellipse_meta_copies(self):
        result = self.reg.to_sky(self.wcs, as_ellipse=True)
        result.meta['text'] = 'new'
        result.visual['color'] = 'green'
        assert result.meta['text'] != self.reg.meta['text']
        assert result.visual['color'] != self.reg.visual['color']


class TestCircleSkyRegionToPixelEllipse:
    """
    Tests for CircleSkyRegion.to_pixel with as_ellipse=True.
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

    def test_to_pixel_as_ellipse(self):
        result = self.reg.to_pixel(self.wcs, as_ellipse=True)
        assert isinstance(result, EllipsePixelRegion)
        assert result.meta == self.meta
        assert result.visual == self.visual

    def test_to_pixel_default_returns_circle(self):
        result = self.reg.to_pixel(self.wcs)
        assert isinstance(result, CirclePixelRegion)

    def test_to_pixel_ellipse_roundtrip(self):
        pix_ellipse = self.reg.to_pixel(self.wcs, as_ellipse=True)
        sky_ellipse = pix_ellipse.to_sky(self.wcs)
        # Roundtrip should recover the original diameter
        assert_quantity_allclose(sky_ellipse.width,
                                 2 * self.radius, rtol=1e-4)
        assert_quantity_allclose(sky_ellipse.height,
                                 2 * self.radius, rtol=1e-4)

    def test_to_pixel_ellipse_center_matches_circle(self):
        pix_circle = self.reg.to_pixel(self.wcs)
        pix_ellipse = self.reg.to_pixel(self.wcs, as_ellipse=True)
        assert_allclose(pix_ellipse.center.x, pix_circle.center.x)
        assert_allclose(pix_ellipse.center.y, pix_circle.center.y)

    def test_to_pixel_as_ellipse_meta_copies(self):
        result = self.reg.to_pixel(self.wcs, as_ellipse=True)
        result.meta['text'] = 'new'
        result.visual['color'] = 'green'
        assert result.meta['text'] != self.reg.meta['text']
        assert result.visual['color'] != self.reg.visual['color']


class TestCircleCovers:
    """
    Test CirclePixelRegion contains and covers boundary semantics using
    explicit expected values.
    """

    @staticmethod
    def setup_circle():
        """
        Create a circle centered at (5, 5) with radius 3.
        """
        return CirclePixelRegion(PixCoord(5, 5), radius=3)

    @pytest.mark.parametrize(('method', 'expected'),
                             [('contains',
                               [True, False, False, False, False]),
                              ('covers',
                               [True, True, False, True, True])])
    def test_array(self, method, expected):
        """
        Test contains/covers boundary semantics for an array mixing
        interior, boundary, and exterior points.
        """
        reg = self.setup_circle()
        # Center (in), right boundary, beyond right (out), top boundary,
        # left boundary
        x = np.array([5, 8, 9, 5, 2])
        y = np.array([5, 5, 5, 8, 5])
        result = getattr(reg, method)(PixCoord(x, y))
        assert_equal(result, np.array(expected))


@pytest.mark.skipif(not HAS_SHAPELY, reason='shapely is required')
class TestCircleShapelyComparison:
    """
    Test that CirclePixelRegion contains and covers match Shapely's
    contains and covers functions (DE-9IM semantics).
    """

    # Test points grouped by location relative to the circle.
    # Shapely approximates a circle with a polygonal ``buffer``, whose
    # vertices land exactly on the cardinal axes. The boundary points
    # below are therefore restricted to the axes (distance 3 from the
    # center) so that they lie on both the true circle and Shapely's
    # polygonal approximation; off-axis boundary points would fall
    # slightly inside the inscribed polygon and disagree.
    boundary_points = [(8, 5), (2, 5), (5, 8), (5, 2)]
    interior_points = [(5, 5), (6, 5), (5, 6), (4, 4)]
    exterior_points = [(9, 5), (1, 5), (5, 9), (5, 1), (0, 0)]
    all_points = boundary_points + interior_points + exterior_points

    @staticmethod
    def setup_circle():
        """
        Create a circle centered at (5, 5) with radius 3.
        """
        return CirclePixelRegion(PixCoord(5, 5), radius=3)

    @staticmethod
    def shapely_circle():
        """
        Create an equivalent Shapely circle (buffer around the center).
        """
        from shapely.geometry import Point as ShapelyPoint
        return ShapelyPoint(5, 5).buffer(3)

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    @pytest.mark.parametrize('point', all_points)
    def test_matches_shapely(self, method, point):
        """
        Test that contains/covers match Shapely for points on the
        boundary, interior, and exterior of the circle.
        """
        from shapely import Point, contains, covers

        shp_func = {'contains': contains, 'covers': covers}[method]
        x, y = point
        reg_result = bool(getattr(self.setup_circle(), method)(PixCoord(x, y)))
        shp_result = bool(shp_func(self.shapely_circle(), Point(x, y)))
        assert reg_result == shp_result
