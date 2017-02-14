# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import pytest, assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from ...core import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB


@pytest.fixture(scope='session')
def wcs():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestCirclePixelRegion:
    def setup(self):
        center = PixCoord(3, 4)
        self.reg = CirclePixelRegion(center, 2)
        self.pixcoord_inside = PixCoord(3, 4)
        self.pixcoord_outside = PixCoord(3, 0)

    def test_repr_str(self):
        reg_repr = '<CirclePixelRegion(PixCoord(x=3, y=4), radius=2)>'
        assert repr(self.reg) == reg_repr

        reg_str = ('Region: CirclePixelRegion\ncenter: PixCoord(x=3, y=4)\n'
                   'radius: 2')
        assert str(self.reg) == reg_str

    def test_contains_scalar(self):
        assert self.reg.contains(self.pixcoord_inside)
        assert self.pixcoord_inside in self.reg

        assert not self.reg.contains(self.pixcoord_outside)
        assert self.pixcoord_outside not in self.reg

    def test_contains_array_1d(self):
        pixcoord = PixCoord([3, 3], [4, 0])
        actual = self.reg.contains(pixcoord)
        expected = [True, False]
        assert_equal(actual, expected)

        with pytest.raises(ValueError) as exc:
            pixcoord in self.reg
        assert 'coord must be scalar' in str(exc)

    def test_contains_array_2d(self):
        pixcoord = PixCoord(
            [[3, 3, 3], [4, 4, 4]],
            [[3, 3, 3], [0, 0, 0]],
        )
        actual = self.reg.contains(pixcoord)
        expected = [[True, True, True], [False, False, False]]
        assert_equal(actual, expected)

    def test_to_mask(self):
        mask = self.reg.to_mask(mode='exact')
        assert_allclose(np.sum(mask.data), 4 * np.pi)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.radius, 2)


class TestCircleSkyRegion:
    def setup(self):
        center = SkyCoord(3 * u.deg, 4 * u.deg)
        self.reg = CircleSkyRegion(center, 2 * u.arcsec)

    def test_repr_str(self):
        if ASTROPY_LT_13:
            reg_repr = ('<CircleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                        'deg\n    (3.0, 4.0)>, radius=2.0 arcsec)>')
            reg_str = ('Region: CircleSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    (3.0, 4.0)>\nradius: 2.0 '
                       'arcsec')
        else:
            reg_repr = ('<CircleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                        'deg\n    ( 3.,  4.)>, radius=2.0 arcsec)>')
            reg_str = ('Region: CircleSkyRegion\ncenter: <SkyCoord (ICRS): '
                       '(ra, dec) in deg\n    ( 3.,  4.)>\nradius: 2.0 '
                       'arcsec')

        assert repr(self.reg) == reg_repr
        assert str(self.reg) == reg_str

    def test_transformation(self, wcs):
        skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
        skycircle = CircleSkyRegion(skycoord, 2 * u.arcsec)

        pixcircle = skycircle.to_pixel(wcs)

        assert_allclose(pixcircle.center.x, -50.5)
        assert_allclose(pixcircle.center.y, 299.5)
        assert_allclose(pixcircle.radius, 0.027777777777828305)

        skycircle2 = pixcircle.to_sky(wcs)

        assert_quantity_allclose(skycircle.center.data.lon, skycircle2.center.data.lon)
        assert_quantity_allclose(skycircle.center.data.lat, skycircle2.center.data.lat)
        assert_quantity_allclose(skycircle2.radius, skycircle.radius)
