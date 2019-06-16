# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS

from ...core import PixCoord
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from ...tests.helpers import make_simple_wcs
from .utils import HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


@pytest.fixture(scope='session')
def wcs():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)

def test_corners():

    xc,yc = 2,2
    angle = 30*u.deg
    width = 2
    height = 1
    reg = RectanglePixelRegion(center=PixCoord(xc, yc),
                               width=width, height=height, angle=angle)

    y1 = yc + np.cos(angle) * height/2 + np.sin(angle) * width/2
    x1 = xc + np.cos(angle) * width/2 - np.sin(angle) * height/2

    assert (x1, y1) in reg.corners

    reg = RectanglePixelRegion(center=PixCoord(xc, yc),
                               width=width, height=height, angle=90*u.deg)
    # simple case: rotate by 90
    np.testing.assert_array_equal([(2.5, 1.), (2.5, 3.), (1.5, 3.), (1.5,1.)],
                                  reg.corners)

    reg = RectanglePixelRegion(center=PixCoord(xc, yc),
                               width=width, height=height, angle=0*u.deg)
    # simpler case: rotate by 0
    np.testing.assert_array_equal([(1, 1.5), (3, 1.5), (3, 2.5), (1, 2.5)],
                                  reg.corners)

    poly = reg.to_polygon()
    assert len(poly.vertices) == 4



class TestRectanglePixelRegion(BaseTestPixelRegion):

    reg = RectanglePixelRegion(center=PixCoord(3, 4), width=4, height=3, angle=5 * u.deg)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 2.5)]
    expected_area = 12
    expected_repr = '<RectanglePixelRegion(PixCoord(x=3, y=4), width=4, height=3, angle=5.0 deg)>'
    expected_str = ('Region: RectanglePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
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
        # Note: `reg.center` is the center, `patch.xy` is the lower-left corner
        assert_allclose(patch.xy, (1.138344, 2.331396), atol=1e-3)

        assert_allclose(patch.get_width(), 4)
        assert_allclose(patch.get_height(), 3)
        # `matplotlib.patches.Rectangle` currently doesn't expose `angle`.
        # See https://github.com/matplotlib/matplotlib/issues/7536
        # In the far future, when it's available in the matplotlib versions
        # we support, we could re-activate a test here.
        # For now, we could also add an assert on `patch.get_verts()` if
        # it's considered important to test that the rotation was done
        # correctly.
        # assert_allclose(patch._angle, 5)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value("deg"), 95)


def test_rectangular_pixel_region_bbox():
    # odd sizes
    width = 7
    height = 3
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0.*u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0.*u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=90.*u.deg)
    assert a.bounding_box.shape == (width, height)

    # even sizes
    width = 8
    height = 4
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0.*u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0.*u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=90.*u.deg)
    assert a.bounding_box.shape == (width, height)


class TestRectangleSkyRegion(BaseTestSkyRegion):

    reg = RectangleSkyRegion(
        center=SkyCoord(3, 4, unit='deg'),
        width=4 * u.deg,
        height=3 * u.deg,
        angle=5 * u.deg,
    )

    expected_repr = ('<RectangleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                'deg\n    ( 3.,  4.)>, width=4.0 deg, height=3.0 '
                'deg, angle=5.0 deg)>')
    expected_str = ('Region: RectangleSkyRegion\ncenter: <SkyCoord '
               '(ICRS): (ra, dec) in deg\n    ( 3.,  4.)>\nwidth: '
               '4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.width.to_value("deg"),  4)
        assert_allclose(reg.height.to_value("deg"), 3)
        assert_allclose(reg.angle.to_value("deg"), 5)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs) == np.array([False, True], dtype='bool'))
