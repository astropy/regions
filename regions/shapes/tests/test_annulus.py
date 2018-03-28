from __future__ import print_function

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..annulus import CircleAnnulusPixelRegion, CircleAnnulusSkyRegion
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion
from .utils import ASTROPY_LT_13


class TestCircleAnnulusPixelRegion(BaseTestPixelRegion):
    reg = CircleAnnulusPixelRegion(PixCoord(3, 4), inner_radius=2, outer_radius=3)
    sample_box = [0, 6, 1, 7]
    inside = [(3, 2)]
    outside = [(3, 0)]
    expected_area = 5 * np.pi
    expected_repr = '<CircleAnnulusPixelRegion(PixCoord(x=3, y=4), inner radius=2, outer radius=3)>'
    expected_str = 'Region: CircleAnnulusPixelRegion\ncenter: PixCoord(x=3, y=4)\ninner radius: 2\nouter radius: 3'

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_radius, 2)
        assert_quantity_allclose(self.reg.outer_radius, 3)

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, CircleAnnulusSkyRegion)


class TestCircleAnnulusSkyRegion(BaseTestSkyRegion):

    reg = CircleAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 20 * u.arcsec, 30 * u.arcsec)
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    if ASTROPY_LT_13:
        expected_repr = ('<CircleAnnulusSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                         'deg\n    (3.0, 4.0)>, inner radius=20.0 arcsec, outer radius=30.0 arcsec)>')
        expected_str = ('Region: CircleAnnulusSkyRegion\ncenter: <SkyCoord (ICRS): '
                        '(ra, dec) in deg\n    (3.0, 4.0)>\ninner radius: 20.0 '
                        'arcsec\nouter radius: 30.0 arcsec')
    else:
        expected_repr = ('<CircleAnnulusSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                         'deg\n    ( 3.,  4.)>, inner radius=20.0 arcsec, outer radius=30.0 arcsec)>')
        expected_str = ('Region: CircleAnnulusSkyRegion\ncenter: <SkyCoord (ICRS): '
                        '(ra, dec) in deg\n    ( 3.,  4.)>\ninner radius: 20.0 '
                        'arcsec\nouter radius: 30.0 arcsec')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_radius, 20*u.arcsec)
        assert_quantity_allclose(self.reg.outer_radius, 30*u.arcsec)

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, CircleAnnulusPixelRegion)
