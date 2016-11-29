# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from numpy.testing import assert_allclose
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import pytest, assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io.fits import getheader
from astropy.wcs import WCS
from ...core import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion
from .utils import ASTROPY_LT_13

try:
    import matplotlib

    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False


class TestCirclePixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = CirclePixelRegion(center, 2)

    def test_str(self):
        expected = 'CirclePixelRegion\ncenter: PixCoord(x=3, y=4)\nradius: 2'
        assert str(self.reg) == expected

    def test_to_mask(self):
        mask = self.reg.to_mask(mode='exact')
        assert_allclose(np.sum(mask.data), 4 * np.pi)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.radius, 2)


class TestCircleSkyRegion:
    def test_str(self):
        center = SkyCoord(3 * u.deg, 4 * u.deg)
        reg = CircleSkyRegion(center, 2 * u.arcsec)

        if ASTROPY_LT_13:
            assert str(
                reg) == 'CircleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n    (3.0, 4.0)>\nradius: 2.0 arcsec'
        else:
            assert str(
                reg) == 'CircleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n    ( 3.,  4.)>\nradius: 2.0 arcsec'

    def test_transformation(self):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        skycircle = CircleSkyRegion(skycoord, 2 * u.arcsec)

        headerfile = get_pkg_data_filename('data/example_header.fits')
        h = getheader(headerfile)
        wcs = WCS(h)

        pixcircle = skycircle.to_pixel(wcs)

        assert_allclose(pixcircle.center.x, -50.5)
        assert_allclose(pixcircle.center.y, 299.5)
        assert_allclose(pixcircle.radius, 0.027777777777828305)

        skycircle2 = pixcircle.to_sky(wcs)

        assert_quantity_allclose(skycircle.center.data.lon, skycircle2.center.data.lon)
        assert_quantity_allclose(skycircle.center.data.lat, skycircle2.center.data.lat)
        assert_quantity_allclose(skycircle2.radius, skycircle.radius)
