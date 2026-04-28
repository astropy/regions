# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
from astropy.coordinates import SkyCoord
from numpy.testing import assert_allclose

from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.core.compound import CompoundPixelRegion, CompoundSkyRegion
from regions.shapes.annulus import (CircleAnnulusPixelRegion,
                                    CircleAnnulusSkyRegion,
                                    EllipseAnnulusPixelRegion,
                                    EllipseAnnulusSkyRegion,
                                    RectangleAnnulusPixelRegion,
                                    RectangleAnnulusSkyRegion)
from regions.shapes.circle import CirclePixelRegion, CircleSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.polygon import (PolygonPixelRegion, PolygonSkyRegion,
                                    RegularPolygonPixelRegion)
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion
from regions.tests.helpers import make_simple_wcs


class TestCirclePixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CirclePixelRegion(PixCoord(3, 4), radius=2, meta=meta,
                            visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.vertices.x) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(n_points=50)
        assert len(poly.vertices.x) == 50

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual
        # meta and visual should be copies
        poly.meta['text'] = 'changed'
        poly.visual['color'] = 'red'
        assert poly.meta['text'] != self.meta['text']
        assert poly.visual['color'] != self.visual['color']

    def test_to_polygon_vertices(self):
        poly = self.reg.to_polygon(n_points=4)
        # theta = [0, pi/2, pi, 3*pi/2]
        assert_allclose(poly.vertices.x, [5, 3, 1, 3], atol=1e-10)
        assert_allclose(poly.vertices.y, [4, 6, 4, 2], atol=1e-10)

    def test_to_polygon_area(self):
        poly = self.reg.to_polygon(n_points=1000)
        assert_allclose(poly.area, self.reg.area, rtol=1e-4)

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        assert poly.contains(PixCoord(3, 4))
        assert not poly.contains(PixCoord(10, 10))


class TestCircleSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = CircleSkyRegion(skycoord, 2 * u.arcsec, meta=meta, visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, PolygonSkyRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon(self.wcs)
        assert len(poly.vertices.ra) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(self.wcs, n_points=50)
        assert len(poly.vertices.ra) == 50


class TestEllipsePixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = EllipsePixelRegion(PixCoord(3, 4), width=6, height=4,
                             angle=0 * u.deg, meta=meta, visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.vertices.x) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(n_points=50)
        assert len(poly.vertices.x) == 50

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual
        poly.meta['text'] = 'changed'
        poly.visual['color'] = 'red'
        assert poly.meta['text'] != self.meta['text']
        assert poly.visual['color'] != self.visual['color']

    def test_to_polygon_vertices_no_rotation(self):
        poly = self.reg.to_polygon(n_points=4)
        # Ellipse: width=6, height=4, center=(3,4), angle=0
        # theta = [0, pi/2, pi, 3*pi/2]
        # x = 0.5*6*cos(theta) + 3 = [6, 3, 0, 3]
        # y = 0.5*4*sin(theta) + 4 = [4, 6, 4, 2]
        assert_allclose(poly.vertices.x, [6, 3, 0, 3], atol=1e-10)
        assert_allclose(poly.vertices.y, [4, 6, 4, 2], atol=1e-10)

    def test_to_polygon_vertices_rotated(self):
        reg = EllipsePixelRegion(PixCoord(0, 0), width=4, height=2,
                                 angle=90 * u.deg)
        poly = reg.to_polygon(n_points=4)
        # After 90 deg rotation, width along y and height along x
        assert_allclose(poly.vertices.x, [0, -1, 0, 1], atol=1e-10)
        assert_allclose(poly.vertices.y, [2, 0, -2, 0], atol=1e-10)

    def test_to_polygon_area(self):
        poly = self.reg.to_polygon(n_points=1000)
        assert_allclose(poly.area, self.reg.area, rtol=1e-4)

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        assert poly.contains(PixCoord(3, 4))
        assert not poly.contains(PixCoord(20, 20))


class TestEllipseSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = EllipseSkyRegion(skycoord, 20 * u.arcsec, 10 * u.arcsec,
                           meta=meta, visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, PolygonSkyRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon(self.wcs)
        assert len(poly.vertices.ra) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(self.wcs, n_points=50)
        assert len(poly.vertices.ra) == 50


class TestRectanglePixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectanglePixelRegion(PixCoord(3, 4), width=4, height=2,
                               angle=0 * u.deg, meta=meta, visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.vertices.x) == 4

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual
        poly.meta['text'] = 'changed'
        poly.visual['color'] = 'red'
        assert poly.meta['text'] != self.meta['text']
        assert poly.visual['color'] != self.visual['color']

    def test_to_polygon_vertices_no_rotation(self):
        poly = self.reg.to_polygon()
        # Rectangle with center=(3,4), width=4, height=2, angle=0
        # Corners: (-2,-1), (2,-1), (2,1), (-2,1) + center
        assert_allclose(sorted(poly.vertices.x), [1, 1, 5, 5])
        assert_allclose(sorted(poly.vertices.y), [3, 3, 5, 5])

    def test_to_polygon_area(self):
        poly = self.reg.to_polygon()
        assert_allclose(poly.area, self.reg.area)

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        assert poly.contains(PixCoord(3, 4))
        assert not poly.contains(PixCoord(20, 20))


class TestRectangleSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = RectangleSkyRegion(skycoord, 20 * u.arcsec, 10 * u.arcsec,
                             meta=meta, visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, PolygonSkyRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon(self.wcs)
        assert len(poly.vertices.ra) == 4


class TestRegularPolygonPixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RegularPolygonPixelRegion(PixCoord(50, 50), nvertices=8,
                                    radius=20, angle=25 * u.deg,
                                    meta=meta, visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.vertices.x) == 8

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual

    def test_to_polygon_vertices(self):
        poly = self.reg.to_polygon()
        assert poly.vertices == self.reg.vertices
        assert poly.origin == PixCoord(0, 0)


class TestCircleAnnulusPixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CircleAnnulusPixelRegion(PixCoord(3, 4), inner_radius=2,
                                   outer_radius=5, meta=meta, visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, CompoundPixelRegion)
        assert isinstance(poly.region1, PolygonPixelRegion)
        assert isinstance(poly.region2, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.region1.vertices.x) == 100
        assert len(poly.region2.vertices.x) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(n_points=50)
        assert len(poly.region1.vertices.x) == 50
        assert len(poly.region2.vertices.x) == 50

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        # Center should be outside (it's an annulus)
        assert not poly.contains(PixCoord(3, 4))
        # Point at radius 3.5 should be inside
        assert poly.contains(PixCoord(6.5, 4))
        # Point far away should be outside
        assert not poly.contains(PixCoord(20, 20))


class TestCircleAnnulusSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = CircleAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                 meta=meta, visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, CompoundSkyRegion)

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(self.wcs, n_points=50)
        assert isinstance(poly, CompoundSkyRegion)


class TestEllipseAnnulusPixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = EllipseAnnulusPixelRegion(PixCoord(3, 4), inner_width=2,
                                    outer_width=6, inner_height=4,
                                    outer_height=8, meta=meta, visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, CompoundPixelRegion)
        assert isinstance(poly.region1, PolygonPixelRegion)
        assert isinstance(poly.region2, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.region1.vertices.x) == 100
        assert len(poly.region2.vertices.x) == 100

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(n_points=50)
        assert len(poly.region1.vertices.x) == 50
        assert len(poly.region2.vertices.x) == 50

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        # Center should be outside (it's an annulus)
        assert not poly.contains(PixCoord(3, 4))
        # Point far away should be outside
        assert not poly.contains(PixCoord(20, 20))


class TestEllipseAnnulusSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = EllipseAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                  50 * u.arcsec, 80 * u.arcsec, meta=meta,
                                  visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, CompoundSkyRegion)

    def test_to_polygon_n_points(self):
        poly = self.reg.to_polygon(self.wcs, n_points=50)
        assert isinstance(poly, CompoundSkyRegion)


class TestRectangleAnnulusPixelRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectangleAnnulusPixelRegion(PixCoord(3, 4), inner_width=2,
                                      outer_width=6, inner_height=4,
                                      outer_height=8, meta=meta,
                                      visual=visual)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon()
        assert isinstance(poly, CompoundPixelRegion)
        assert isinstance(poly.region1, PolygonPixelRegion)
        assert isinstance(poly.region2, PolygonPixelRegion)

    def test_to_polygon_nvertices(self):
        poly = self.reg.to_polygon()
        assert len(poly.region1.vertices.x) == 4
        assert len(poly.region2.vertices.x) == 4

    def test_to_polygon_meta(self):
        poly = self.reg.to_polygon()
        assert poly.meta == self.meta
        assert poly.visual == self.visual

    def test_to_polygon_contains(self):
        poly = self.reg.to_polygon()
        # Center should be outside (it's an annulus)
        assert not poly.contains(PixCoord(3, 4))
        # Point far away should be outside
        assert not poly.contains(PixCoord(20, 20))


class TestRectangleAnnulusSkyRegionToPolygon:
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = RectangleAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                    50 * u.arcsec, 80 * u.arcsec, meta=meta,
                                    visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_to_polygon_type(self):
        poly = self.reg.to_polygon(self.wcs)
        assert isinstance(poly, CompoundSkyRegion)
