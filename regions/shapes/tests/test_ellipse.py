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
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestEllipsePixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = EllipsePixelRegion(center=PixCoord(3, 4), width=4, height=3,
                             angle=5 * u.deg, meta=meta, visual=visual)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 4)]
    expected_area = 3 * np.pi
    expected_repr = ('<EllipsePixelRegion(center=PixCoord(x=3, y=4), '
                     'width=4, height=3, angle=5.0 deg)>')
    expected_str = ('Region: EllipsePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
                    'width: 4\nheight: 3\nangle: 5.0 deg')

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
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.width, 4)
        assert_allclose(patch.height, 3)
        assert_allclose(patch.angle, 5)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value('deg'), 95)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.width = 3
        assert reg != self.reg

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

    def test_region_bbox_zero_size(self):
        with pytest.raises(ValueError):
            EllipsePixelRegion(PixCoord(50, 50), width=0, height=0,
                               angle=0. * u.deg)

        with pytest.raises(ValueError):
            EllipsePixelRegion(PixCoord(50, 50), width=10, height=0,
                               angle=0. * u.deg)

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    @pytest.mark.parametrize('sync', (False, True))
    def test_as_mpl_selector(self, sync):
        import matplotlib.pyplot as plt

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        _, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(
                mode='subpixels', subpixels=10).to_image(data.shape)

        # For now this will only work with unrotated ellipses. Once this
        # works with rotated ellipses, the following exception check can
        # be removed as well as the ``angle=0 * u.deg`` in the call to
        # copy() below.
        with pytest.raises(NotImplementedError,
                           match=('Cannot create matplotlib selector'
                                  ' for rotated ellipse.')):
            self.reg.as_mpl_selector(ax)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask, sync=sync)  # noqa: E501, F841

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

            assert_equal(mask, region.to_mask(
                mode='subpixels', subpixels=10).to_image(data.shape))

        else:

            assert_allclose(region.center.x, 3)
            assert_allclose(region.center.y, 4)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, 0)

        match = 'Cannot attach more than one selector to a reg'
        with pytest.raises(AttributeError, match=match):
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

        _, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(
                mode='subpixels', subpixels=10).to_image(data.shape)

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

        # For drag_from_anywhere=False this will have created a new 1x1
        # rectangle.
        if anywhere:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
        else:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)

        assert_equal(mask, region.to_mask(
            mode='subpixels', subpixels=10).to_image(data.shape))

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

        _, ax = plt.subplots()
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(
                mode='subpixels', subpixels=10).to_image(data.shape)

        region = self.reg.copy(angle=0 * u.deg)
        region.visual = {'color': 'red'}

        if 'twit' in userargs:
            match = r'__init__.. got an unexpected keyword argument'
            with pytest.raises(TypeError, match=match):
                selector = region.as_mpl_selector(
                    ax, callback=update_mask, **userargs)
        else:
            selector = region.as_mpl_selector(
                ax, callback=update_mask, **userargs)
            assert region._mpl_selector.artists[0].get_edgecolor() == (
                1, 0, 0, 1)

            if 'props' in userargs:
                assert region._mpl_selector.artists[0].get_facecolor() == (
                    0, 0, 1, 1)
                assert region._mpl_selector.artists[0].get_linewidth() == 2
            else:
                assert region._mpl_selector.artists[0].get_facecolor() == (
                    0, 0, 0, 0)
                assert region._mpl_selector.artists[0].get_linewidth() == 1

                for key, val in userargs.items():
                    assert getattr(region._mpl_selector, key) == val
                    assert getattr(selector, key) == val

        plt.close()


class TestEllipseSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = EllipseSkyRegion(center=SkyCoord(3, 4, unit='deg'), width=4 * u.deg,
                           height=3 * u.deg, angle=5 * u.deg, meta=meta,
                           visual=visual)

    expected_repr = ('<EllipseSkyRegion(center=<SkyCoord (ICRS): (ra, dec) '
                     'in deg\n    (3., 4.)>, width=4.0 deg, height=3.0 deg,'
                     ' angle=5.0 deg)>')
    expected_str = ('Region: EllipseSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>\nwidth: 4.0 deg\n'
                    'height: 3.0 deg\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.width.to_value('deg'), 4)
        assert_allclose(reg.height.to_value('deg'), 3)
        assert_allclose(reg.angle.to_value('deg'), 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_dimension_center(self):
        center = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        width = 2 * u.arcsec
        height = 3 * u.arcsec
        with pytest.raises(ValueError) as excinfo:
            EllipseSkyRegion(center, width, height)
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
        reg.width = 3 * u.deg
        assert reg != self.reg


class TestEllipseCovers:
    """
    Test EllipsePixelRegion contains and covers boundary semantics using
    explicit expected values.
    """

    @staticmethod
    def setup_ellipse():
        """
        Create an axis-aligned ellipse centered at (5, 5).
        """
        # width=6 (semi-major=3 along x), height=4 (semi-minor=2 along y)
        return EllipsePixelRegion(PixCoord(5, 5), width=6, height=4,
                                  angle=0 * u.deg)

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
        reg = self.setup_ellipse()
        # Center (in), right boundary, beyond right (out), top boundary,
        # left boundary
        x = np.array([5, 8, 9, 5, 2])
        y = np.array([5, 5, 5, 7, 5])
        result = getattr(reg, method)(PixCoord(x, y))
        assert_equal(result, np.array(expected))

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    def test_rotated_ellipse(self, method):
        """
        Test boundary semantics for an ellipse rotated by 90 degrees,
        which swaps the effect of the semi-major and semi-minor axes.
        """
        # After 90 deg rotation, the semi-major axis (3) is vertical and
        # the semi-minor axis (2) is horizontal, so the boundary points
        # are (5, 8), (5, 2) (top/bottom) and (7, 5), (3, 5) (left/right).
        reg = EllipsePixelRegion(PixCoord(5, 5), width=6, height=4,
                                 angle=90 * u.deg)
        boundary = [(5, 8), (5, 2), (7, 5), (3, 5)]
        # Boundary points are excluded by contains but included by covers.
        on_boundary = (method == 'covers')

        for x, y in boundary:
            assert bool(getattr(reg, method)(PixCoord(x, y))) is on_boundary

        # Interior points are always inside; exterior points never are.
        assert getattr(reg, method)(PixCoord(5, 5))  # center
        assert getattr(reg, method)(PixCoord(5, 7))  # inside vertically
        assert not getattr(reg, method)(PixCoord(5, 9))  # beyond top
        assert not getattr(reg, method)(PixCoord(8, 5))  # beyond right


@pytest.mark.skipif(not HAS_SHAPELY, reason='shapely is required')
class TestEllipseShapelyComparison:
    """
    Test that EllipsePixelRegion contains and covers match Shapely's
    contains and covers functions (DE-9IM semantics).
    """

    # Test points grouped by location relative to the ellipse.
    # Shapely approximates the ellipse with a scaled polygonal
    # ``buffer``, whose vertices land on the semi-major/semi-minor axes.
    # The boundary points below are therefore restricted to those axes
    # so that they lie on both the true ellipse and Shapely's polygonal
    # approximation; off-axis boundary points would fall slightly inside
    # the inscribed polygon and disagree.
    boundary_points = [(8, 5), (2, 5), (5, 7), (5, 3)]
    interior_points = [(5, 5), (6, 5), (5, 6), (4, 4), (7, 5)]
    exterior_points = [(9, 5), (1, 5), (5, 8), (5, 2), (0, 0), (8, 6)]
    all_points = boundary_points + interior_points + exterior_points

    @staticmethod
    def setup_ellipse():
        """
        Create an axis-aligned ellipse centered at (5, 5).
        """
        return EllipsePixelRegion(PixCoord(5, 5), width=6, height=4,
                                  angle=0 * u.deg)

    @staticmethod
    def shapely_ellipse():
        """
        Create an equivalent Shapely ellipse by scaling a unit circle.
        """
        from shapely import affinity
        from shapely.geometry import Point as ShapelyPoint

        # width=6 -> semi-major a=3, height=4 -> semi-minor b=2
        circle = ShapelyPoint(5, 5).buffer(1)
        return affinity.scale(circle, xfact=3, yfact=2)

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    @pytest.mark.parametrize('point', all_points)
    def test_matches_shapely(self, method, point):
        """
        Test that contains/covers match Shapely for points on the
        boundary, interior, and exterior of the ellipse.
        """
        from shapely import Point, contains, covers

        shp_func = {'contains': contains, 'covers': covers}[method]
        x, y = point
        reg_result = bool(
            getattr(self.setup_ellipse(), method)(PixCoord(x, y)))
        shp_result = bool(shp_func(self.shapely_ellipse(), Point(x, y)))
        assert reg_result == shp_result
