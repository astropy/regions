# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion
from .utils import HAS_MATPLOTLIB  # noqa


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestEllipsePixelRegion(BaseTestPixelRegion):
    reg = EllipsePixelRegion(center=PixCoord(3, 4), width=4, height=3,
                             angle=5 * u.deg)
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

    def test_region_bbox_zero_size(self):
        reg = EllipsePixelRegion(PixCoord(50, 50), width=0, height=0,
                                 angle=0. * u.deg)
        assert reg.bounding_box.shape == (1, 1)

        reg = EllipsePixelRegion(PixCoord(50, 50), width=10, height=0,
                                 angle=0. * u.deg)
        assert reg.bounding_box.shape == (1, 11)

        reg = EllipsePixelRegion(PixCoord(50, 50), width=0, height=10,
                                 angle=0. * u.deg)
        assert reg.bounding_box.shape == (11, 1)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value("deg"), 95)

    # TODO: Is this MatplotlibDeprecationWarning something to worry about?
    @pytest.mark.filterwarnings(r"ignore:The 'rectprops' parameter of "
                                r"__init__\(\) has been renamed 'props'")
    @pytest.mark.parametrize('sync', (False, True))
    def test_as_mpl_selector(self, sync):

        plt = pytest.importorskip('matplotlib.pyplot')

        data = np.random.random((16, 16))
        mask = np.zeros_like(data)

        ax = plt.subplot(1, 1, 1)
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels',
                                  subpixels=10).to_image(data.shape)

        # For now this will only work with unrotated ellipses. Once this
        # works with rotated ellipses, the following exception check can
        # be removed as well as the ``angle=0 * u.deg`` in the call to
        # copy() below.
        with pytest.raises(NotImplementedError,
                           match=('Cannot create matplotlib selector for '
                                  'rotated ellipse.')):
            self.reg.as_mpl_selector(ax)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask, sync=sync)  # noqa

        from matplotlib.backend_bases import MouseEvent, MouseButton

        x, y = ax.transData.transform([[7.3, 4.4]])[0]
        ax.figure.canvas.callbacks.process('button_press_event',
                                           MouseEvent('button_press_event',
                                                      ax.figure.canvas, x, y,
                                                      button=MouseButton.LEFT))
        x, y = ax.transData.transform([[9.3, 5.4]])[0]
        ax.figure.canvas.callbacks.process('motion_notify_event',
                                           MouseEvent('button_press_event',
                                                      ax.figure.canvas, x, y,
                                                      button=MouseButton.LEFT))
        x, y = ax.transData.transform([[9.3, 5.4]])[0]
        ax.figure.canvas.callbacks.process('button_release_event',
                                           MouseEvent('button_press_event',
                                                      ax.figure.canvas, x, y,
                                                      button=MouseButton.LEFT))

        ax.figure.canvas.draw()

        if sync:

            assert_allclose(region.center.x, 8.3)
            assert_allclose(region.center.y, 4.9)
            assert_allclose(region.width, 2)
            assert_allclose(region.height, 1)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask,
                         region.to_mask(mode='subpixels',
                                        subpixels=10).to_image(data.shape))

        else:

            assert_allclose(region.center.x, 3)
            assert_allclose(region.center.y, 4)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, 0)

        with pytest.raises(Exception, match=('Cannot attach more than one '
                                             'selector to a region.')):
            region.as_mpl_selector(ax)


class TestEllipseSkyRegion(BaseTestSkyRegion):
    reg = EllipseSkyRegion(
        center=SkyCoord(3, 4, unit='deg'),
        width=4 * u.deg,
        height=3 * u.deg,
        angle=5 * u.deg,
    )

    expected_repr = ('<EllipseSkyRegion(center=<SkyCoord (ICRS): (ra, dec) '
                     'in deg\n    (3., 4.)>, width=4.0 deg, height=3.0 deg,'
                     ' angle=5.0 deg)>')
    expected_str = ('Region: EllipseSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>\nwidth: 4.0 deg\n'
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
        with pytest.raises(ValueError) as excinfo:
            EllipseSkyRegion(center, width, height)
        estr = 'The center must be a scalar SkyCoord object'
        assert estr in str(excinfo.value)

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, True], dtype='bool'))
