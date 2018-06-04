# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..circle import CirclePixelRegion, CircleSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


@pytest.fixture(scope='session')
def wcs():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestCirclePixelRegion(BaseTestPixelRegion):

    reg = CirclePixelRegion(PixCoord(3, 4), radius=2)
    sample_box = [0, 6, 1, 7]
    inside = [(3, 4)]
    outside = [(3, 0)]
    expected_area = 4 * np.pi
    expected_repr = '<CirclePixelRegion(PixCoord(x=3, y=4), radius=2)>'
    expected_str = 'Region: CirclePixelRegion\ncenter: PixCoord(x=3, y=4)\nradius: 2'

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.radius, self.reg.radius)

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        assert_allclose(patch.center, (3, 4))
        assert_allclose(patch.radius, 2)


class TestCircleSkyRegion(BaseTestSkyRegion):

    reg = CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 2 * u.arcsec)

    if ASTROPY_LT_13:
        expected_repr = ('<CircleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                         'deg\n    (3.0, 4.0)>, radius=2.0 arcsec)>')
        expected_str = ('Region: CircleSkyRegion\ncenter: <SkyCoord (ICRS): '
                        '(ra, dec) in deg\n    (3.0, 4.0)>\nradius: 2.0 '
                        'arcsec')
    else:
        expected_repr = ('<CircleSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                         'deg\n    ( 3.,  4.)>, radius=2.0 arcsec)>')
        expected_str = ('Region: CircleSkyRegion\ncenter: <SkyCoord (ICRS): '
                        '(ra, dec) in deg\n    ( 3.,  4.)>\nradius: 2.0 '
                        'arcsec')

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

    def test_dimension_center(self):
        center = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        radius = 2 * u.arcsec
        with pytest.raises(ValueError) as err:
            CircleSkyRegion(center, radius)
        assert 'the centre should be a 0D SkyCoord object' in str(err)
