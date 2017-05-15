# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

from numpy.testing import assert_allclose

import astropy.units as u
from astropy.tests.helper import pytest
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import SkyCoord

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB  # noqa


class TestEllipsePixelRegion():

    def setup(self):
        center = PixCoord(3, 4)
        self.reg = EllipsePixelRegion(
            center=center,
            width=4,
            height=3,
            angle=5 * u.deg,
        )

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.width, self.reg.width)
        assert_allclose(reg_new.height, self.reg.height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)

    def test_repr_str(self):
        reg_repr = ('<EllipsePixelRegion(PixCoord(x=3, y=4), width=4, height=3'
                    ', angle=5.0 deg)>')
        assert repr(self.reg) == reg_repr

        reg_str = ('Region: EllipsePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
                   'width: 4\nheight: 3\nangle: 5.0 deg')
        assert str(self.reg) == reg_str

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.width, 8)
        assert_allclose(patch.height, 6)
        assert_allclose(patch.angle, 5)


class TestEllipseSkyRegion:
    def setup(self):
        center = SkyCoord(3, 4, unit='deg')
        self.reg = EllipseSkyRegion(
            center=center,
            width=4 * u.deg,
            height=3 * u.deg,
            angle=5 * u.deg,
        )

    def test_repr_str(self):
        if ASTROPY_LT_13:
            reg_repr = ('<EllipseSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                        'deg\n    (3.0, 4.0)>, width=4.0 deg, height=3.0 deg,'
                        ' angle=5.0 deg)>')
            reg_str = ('Region: EllipseSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    (3.0, 4.0)>\nwidth: 4.0 deg\n'
                       'height: 3.0 deg\nangle: 5.0 deg')
        else:
            reg_repr = ('<EllipseSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                        'deg\n    ( 3.,  4.)>, width=4.0 deg, height=3.0 deg,'
                        ' angle=5.0 deg)>')
            reg_str = ('Region: EllipseSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    ( 3.,  4.)>\nwidth: 4.0 deg\n'
                       'height: 3.0 deg\nangle: 5.0 deg')

        assert repr(self.reg) == reg_repr
        assert str(self.reg) == reg_str


def test_ellipse_pixel_region_bbox():
    a = 7
    b = 3
    reg = EllipsePixelRegion(PixCoord(50, 50), width=a, height=b,
                             angle=0.*u.deg)
    assert reg.bounding_box.shape == (2*b + 1, 2*a + 1)

    reg = EllipsePixelRegion(PixCoord(50.5, 50.5), width=a, height=b,
                             angle=0.*u.deg)
    assert reg.bounding_box.shape == (2*b, 2*a)

    reg = EllipsePixelRegion(PixCoord(50, 50), width=a, height=b,
                             angle=90.*u.deg)
    assert reg.bounding_box.shape == (2*a + 1, 2*b + 1)
