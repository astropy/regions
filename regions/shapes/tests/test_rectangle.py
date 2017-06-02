# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from numpy.testing import assert_allclose
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from ...core import PixCoord
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from ...tests.helpers import make_simple_wcs
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


class TestRectanglePixelRegion(BaseTestPixelRegion):

    reg = RectanglePixelRegion(center=PixCoord(3, 4), width=4, height=3, angle=5 * u.deg)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 2.5)]
    expected_area = 12
    expected_repr = '<RectanglePixelRegion(PixCoord(x=3, y=4), width=4, height=3, angle=5.0 deg)>'
    expected_str = ('Region: RectanglePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
                    'width: 4\nheight: 3\nangle: 5.0 deg')

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.width, self.reg.width)
        assert_allclose(reg_new.height, self.reg.height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
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

    if ASTROPY_LT_13:
        expected_repr = ('<RectangleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                    'deg\n    (3.0, 4.0)>, width=4.0 deg, height=3.0 '
                    'deg, angle=5.0 deg)>')
        expected_str = ('Region: RectangleSkyRegion\ncenter: <SkyCoord '
                   '(ICRS): (ra, dec) in deg\n    (3.0, 4.0)>\nwidth: '
                   '4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')
    else:
        expected_repr = ('<RectangleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                    'deg\n    ( 3.,  4.)>, width=4.0 deg, height=3.0 '
                    'deg, angle=5.0 deg)>')
        expected_str = ('Region: RectangleSkyRegion\ncenter: <SkyCoord '
                   '(ICRS): (ra, dec) in deg\n    ( 3.,  4.)>\nwidth: '
                   '4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')
