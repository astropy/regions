# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose, assert_equal

from regions._utils.optional_deps import HAS_MATPLOTLIB, HAS_SHAPELY
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


def test_corners():
    xc, yc = 2, 2
    angle = 30 * u.deg
    width = 2
    height = 1
    reg = RectanglePixelRegion(PixCoord(xc, yc), width=width, height=height,
                               angle=angle)

    y1 = yc + np.cos(angle) * (height / 2) + np.sin(angle) * (width / 2)
    x1 = xc + np.cos(angle) * (width / 2) - np.sin(angle) * (height / 2)

    assert (x1, y1) in reg.corners

    reg = RectanglePixelRegion(PixCoord(xc, yc), width=width, height=height,
                               angle=90 * u.deg)
    # simple case: rotate by 90
    np.testing.assert_allclose([(2.5, 1.), (2.5, 3.), (1.5, 3.), (1.5, 1.)],
                               reg.corners)

    reg = RectanglePixelRegion(center=PixCoord(xc, yc), width=width,
                               height=height, angle=0 * u.deg)
    # simpler case: rotate by 0
    np.testing.assert_array_equal([(1, 1.5), (3, 1.5), (3, 2.5), (1, 2.5)],
                                  reg.corners)

    poly = reg.to_polygon()
    assert len(poly.vertices) == 4


class TestRectanglePixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectanglePixelRegion(center=PixCoord(3, 4), width=4, height=3,
                               angle=5 * u.deg, meta=meta, visual=visual)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 2.5)]
    expected_area = 12
    expected_repr = ('<RectanglePixelRegion(center=PixCoord(x=3, y=4), '
                     'width=4, height=3, angle=5.0 deg)>')
    expected_str = ('Region: RectanglePixelRegion\ncenter: PixCoord(x=3, '
                    'y=4)\nwidth: 4\nheight: 3\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.width == 4
        assert reg.height == 3
        assert_allclose(reg.angle.to_value('deg'), 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.width, self.reg.width)
        assert_allclose(reg_new.height, self.reg.height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)
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
        # Note: `reg.center` is the center, `patch.xy` is the lower-left
        # corner
        assert_allclose(patch.xy, (1.138344, 2.331396), atol=1e-3)
        assert_allclose(patch.get_width(), 4)
        assert_allclose(patch.get_height(), 3)
        assert_allclose(patch.angle, 5)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value('deg'), 95)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.angle = 35 * u.deg
        assert reg != self.reg

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    @pytest.mark.parametrize('sync', (False, True))
    def test_as_mpl_selector(self, sync):
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        fig, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        # For now this will only work with unrotated rectangles. Once
        # this works with rotated rectangles, the following exception
        # check can be removed as well as the ``angle=0 * u.deg`` in the
        # call to copy() below.
        with pytest.raises(NotImplementedError,
                           match=('Cannot create matplotlib selector for rotated rectangle.')):
            self.reg.as_mpl_selector(ax)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask, sync=sync)  # noqa: F841

        from matplotlib.backend_bases import MouseEvent
        canvas = ax.figure.canvas

        evt = MouseEvent('button_press_event', canvas,
                         *ax.transData.transform((7.3, 4.4)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('motion_notify_event', canvas,
                         *ax.transData.transform((9.3, 5.4)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('button_release_event', canvas,
                         *ax.transData.transform((9.3, 5.4)), button=1)
        canvas.callbacks.process(evt.name, evt)
        ax.figure.canvas.draw()

        if sync:
            assert_allclose(region.center.x, 8.3)
            assert_allclose(region.center.y, 4.9)
            assert_allclose(region.width, 2)
            assert_allclose(region.height, 1)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, region.to_mask(mode='subpixels', subpixels=10).to_image(data.shape))

        else:
            assert_allclose(region.center.x, 3)
            assert_allclose(region.center.y, 4)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, 0)

        with pytest.raises(AttributeError, match=('Cannot attach more than one selector to a reg')):
            region.as_mpl_selector(ax)

        plt.close()

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    @pytest.mark.parametrize('anywhere', (False, True))
    def test_mpl_selector_drag(self, anywhere):
        """
        Test dragging of entire region from central handle and anywhere.
        """
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        fig, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask,
                                          drag_from_anywhere=anywhere)
        assert selector.drag_from_anywhere is anywhere
        assert region._mpl_selector.drag_from_anywhere is anywhere

        # click_and_drag(selector, start=(3, 4), end=(3.5, 4.5))

        from matplotlib.backend_bases import MouseEvent
        canvas = ax.figure.canvas

        evt = MouseEvent('button_press_event', canvas,
                         *ax.transData.transform((3, 4)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('motion_notify_event', canvas,
                         *ax.transData.transform((3.5, 4.5)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('button_release_event', canvas,
                         *ax.transData.transform((3.5, 4.5)), button=1)
        canvas.callbacks.process(evt.name, evt)
        ax.figure.canvas.draw()

        assert_allclose(region.center.x, 3.5)
        assert_allclose(region.center.y, 4.5)
        assert_allclose(region.width, 4)
        assert_allclose(region.height, 3)

        evt = MouseEvent('button_press_event', canvas,
                         *ax.transData.transform((3.25, 4.25)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('motion_notify_event', canvas,
                         *ax.transData.transform((4.25, 5.25)), button=1)
        canvas.callbacks.process(evt.name, evt)
        evt = MouseEvent('button_release_event', canvas,
                         *ax.transData.transform((4.25, 5.25)), button=1)
        canvas.callbacks.process(evt.name, evt)
        ax.figure.canvas.draw()

        # For drag_from_anywhere=False this will have created a new 1x1 rectangle.
        if anywhere:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
        else:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)

        assert_equal(mask, region.to_mask(mode='subpixels', subpixels=10).to_image(data.shape))

        plt.close()

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    @pytest.mark.parametrize('userargs',
                             ({'useblit': True},
                              {'grab_range': 20, 'minspanx': 5, 'minspany': 4},
                              {'props': {'facecolor': 'blue', 'linewidth': 2}},
                              {'twit': 'gumby'}))
    def test_mpl_selector_kwargs(self, userargs):
        """
        Test that additional kwargs are passed to selector.
        """
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        fig, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        region = self.reg.copy(angle=0 * u.deg)

        if 'twit' in userargs:
            with pytest.raises(TypeError, match=(r'__init__.. got an unexpected keyword argument')):
                selector = region.as_mpl_selector(ax, callback=update_mask, **userargs)
        else:
            selector = region.as_mpl_selector(ax, callback=update_mask, **userargs)
            assert region._mpl_selector.artists[0].get_edgecolor() == (0, 0, 1, 1)

            if 'props' in userargs:
                assert region._mpl_selector.artists[0].get_facecolor() == (0, 0, 1, 1)
                assert region._mpl_selector.artists[0].get_linewidth() == 2
            else:
                assert region._mpl_selector.artists[0].get_facecolor() == (0, 0, 0, 0)
                assert region._mpl_selector.artists[0].get_linewidth() == 1

                for key, val in userargs.items():
                    assert getattr(region._mpl_selector, key) == val
                    assert getattr(selector, key) == val

        plt.close()


def test_rectangular_pixel_region_bbox():
    # odd sizes
    width = 7
    height = 3
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=90. * u.deg)
    assert a.bounding_box.shape == (width, height)

    # even sizes
    width = 8
    height = 4
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=90. * u.deg)
    assert a.bounding_box.shape == (width, height)


class TestRectangleSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectangleSkyRegion(center=SkyCoord(3, 4, unit='deg'),
                             width=4 * u.deg, height=3 * u.deg,
                             angle=5 * u.deg, meta=meta, visual=visual)

    expected_repr = ('<RectangleSkyRegion(center=<SkyCoord (ICRS): (ra, dec) '
                     'in deg\n    (3., 4.)>, width=4.0 deg, height=3.0 deg, '
                     'angle=5.0 deg)>')
    expected_str = ('Region: RectangleSkyRegion\ncenter: <SkyCoord '
                    '(ICRS): (ra, dec) in deg\n    (3., 4.)>\nwidth: '
                    '4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.width.to_value('deg'), 4)
        assert_allclose(reg.height.to_value('deg'), 3)
        assert_allclose(reg.angle.to_value('deg'), 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, True], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.angle = 10 * u.deg
        assert reg != self.reg


class TestRectangleCovers:
    """
    Test RectanglePixelRegion contains and covers boundary semantics
    using explicit expected values.
    """

    @staticmethod
    def setup_rectangle():
        """
        Create an axis-aligned rectangle centered at (5, 5).
        """
        # width=6, height=4, so corners at (2, 3), (8, 3), (8, 7), (2, 7)
        return RectanglePixelRegion(PixCoord(5, 5), width=6, height=4,
                                    angle=0 * u.deg)

    @pytest.mark.parametrize(('method', 'expected'),
                             [('contains', [True, False, False, False, False]),
                              ('covers', [True, True, False, True, True])])
    def test_array(self, method, expected):
        """
        Test contains/covers boundary semantics for an array mixing
        interior, boundary, and exterior points.
        """
        reg = self.setup_rectangle()
        # Center (in), vertex, beyond right (out), bottom edge, right edge
        x = np.array([5, 2, 9, 5, 8])
        y = np.array([5, 3, 5, 3, 5])
        result = getattr(reg, method)(PixCoord(x, y))
        assert_equal(result, np.array(expected))

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    def test_rotated_rectangle(self, method):
        """
        Test boundary semantics for a rectangle rotated by 45 degrees.

        Exact-boundary comparisons are sensitive to floating-point
        rounding, so points are nudged just inside and just outside the
        corner along the center-corner direction.
        """
        reg = RectanglePixelRegion(PixCoord(5, 5), width=4, height=2,
                                   angle=45 * u.deg)

        # Center is always inside; a far point is always outside.
        assert getattr(reg, method)(PixCoord(5, 5))
        assert not getattr(reg, method)(PixCoord(10, 10))

        cos45 = np.cos(np.radians(45))
        sin45 = np.sin(np.radians(45))
        hw, hh = 2, 1  # half width, half height
        cx = 5 + cos45 * hw - sin45 * hh
        cy = 5 + sin45 * hw + cos45 * hh

        inside = PixCoord(5 + 0.99 * (cx - 5), 5 + 0.99 * (cy - 5))
        outside = PixCoord(5 + 1.01 * (cx - 5), 5 + 1.01 * (cy - 5))
        assert getattr(reg, method)(inside)
        assert not getattr(reg, method)(outside)


@pytest.mark.skipif(not HAS_SHAPELY, reason='shapely is required')
class TestRectangleShapelyComparison:
    """
    Test that RectanglePixelRegion contains and covers match Shapely's
    contains and covers functions (DE-9IM semantics).
    """

    # Test points grouped by location relative to the rectangle
    vertices = [(2, 3), (8, 3), (8, 7), (2, 7)]
    edge_points = [(5, 3), (8, 5), (5, 7), (2, 5)]
    interior_points = [(5, 5), (3, 4), (7, 6), (4, 5)]
    exterior_points = [(0, 0), (1, 5), (9, 5), (5, 2), (5, 8)]
    all_points = vertices + edge_points + interior_points + exterior_points

    @staticmethod
    def setup_rectangle():
        """
        Create an axis-aligned rectangle centered at (5, 5).
        """
        return RectanglePixelRegion(PixCoord(5, 5), width=6, height=4,
                                    angle=0 * u.deg)

    @staticmethod
    def shapely_rectangle():
        """
        Create an equivalent Shapely rectangle polygon.
        """
        from shapely import Polygon
        return Polygon([(2, 3), (8, 3), (8, 7), (2, 7)])

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    @pytest.mark.parametrize('point', all_points)
    def test_matches_shapely(self, method, point):
        """
        Test that contains/covers match Shapely for points on the
        vertices, edges, interior, and exterior of the rectangle.
        """
        from shapely import Point, contains, covers

        shp_func = {'contains': contains, 'covers': covers}[method]
        x, y = point
        reg = self.setup_rectangle()
        reg_result = bool(getattr(reg, method)(PixCoord(x, y)))
        shp_result = bool(shp_func(self.shapely_rectangle(), Point(x, y)))
        assert reg_result == shp_result
