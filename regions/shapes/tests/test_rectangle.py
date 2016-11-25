# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from .utils import ASTROPY_LT_13


def test_ellipse_pixel():
    center = PixCoord(3, 4)
    reg = RectanglePixelRegion(center, 3, 4, 5 * u.deg)

    assert str(reg) == 'RectanglePixelRegion\ncenter: PixCoord(x=3, y=4)\nheight: 3\nwidth: 4\nangle: 5.0 deg'


def test_ellipse_sky():
    center = SkyCoord(3, 4, unit='deg')
    reg = RectangleSkyRegion(center, 3 * u.deg, 4 * u.deg, 5 * u.deg)
    if ASTROPY_LT_13:
        expected = ('RectangleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                    '    (3.0, 4.0)>\nheight: 3.0 deg\nwidth: 4.0 deg\nangle: 5.0 deg')
    else:
        expected = ('RectangleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n'
                    '    ( 3.,  4.)>\nheight: 3.0 deg\nwidth: 4.0 deg\nangle: 5.0 deg')
    assert str(reg) == expected
