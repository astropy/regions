from __future__ import print_function

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..annulus import *
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


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

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.inner_radius == 2
        assert reg.outer_radius == 3
        assert reg.visual == {}
        assert reg.meta == {}

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, CircleAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))


class TestCircleAnnulusSkyRegion(BaseTestSkyRegion):
    reg = CircleAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 20 * u.arcsec, 30 * u.arcsec)
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<CircleAnnulusSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    ( 3.,  4.)>, inner radius=20.0 arcsec, outer radius=30.0 arcsec)>')
    expected_str = ('Region: CircleAnnulusSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    ( 3.,  4.)>\ninner radius: 20.0 '
                    'arcsec\nouter radius: 30.0 arcsec')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_radius, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.outer_radius, 30 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_radius.to_value("arcsec"), 20)
        assert_allclose(reg.outer_radius.to_value("arcsec"), 30)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, CircleAnnulusPixelRegion)


class TestEllipseAnnulusPixelRegion(BaseTestPixelRegion):
    reg = EllipseAnnulusPixelRegion(PixCoord(3, 4), 2, 5, 5, 8)
    sample_box = [1, 6, 0, 9]
    inside = [(3, 7)]
    outside = [(3, 4)]
    expected_area = 7.5 * np.pi
    expected_repr = ('<EllipseAnnulusPixelRegion(PixCoord(x=3, y=4), '
                     'inner width=2, inner height=5, outer width=5, '
                     'outer height=8, angle=0.0 deg)>')
    expected_str = ('Region: EllipseAnnulusPixelRegion\ncenter: PixCoord(x=3, y=4)'
                    '\ninner width: 2\ninner height: 5\nouter width: 5\n'
                    'outer height: 8\nangle: 0.0 deg')

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_width, 2)
        assert_quantity_allclose(self.reg.inner_height, 5)
        assert_quantity_allclose(self.reg.outer_width, 5)
        assert_quantity_allclose(self.reg.outer_height, 8)

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.inner_width == 2
        assert reg.inner_height == 5
        assert reg.outer_width == 5
        assert reg.outer_height == 8
        assert_allclose(reg.angle.to_value("deg"), 0)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, EllipseAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value("deg"), 90)


class TestEllipseAnnulusSkyRegion(BaseTestSkyRegion):
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    reg = EllipseAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                  50 * u.arcsec, 80 * u.arcsec)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<EllipseAnnulusSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    (3., 4.)>, inner width=20.0 arcsec, '
                     'inner height=50.0 arcsec, outer width=50.0 arcsec, '
                     'outer height=80.0 arcsec, angle=0.0 deg)>')
    expected_str = ('Region: EllipseAnnulusSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>\ninner width: 20.0 arcsec\n'
                    'inner height: 50.0 arcsec\nouter width: 50.0 arcsec\n'
                    'outer height: 80.0 arcsec\nangle: 0.0 deg')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_width, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.inner_height, 50 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_width.to_value("arcsec"), 20)
        assert_allclose(reg.inner_height.to_value("arcsec"), 50)
        assert_allclose(reg.outer_width.to_value("arcsec"), 50)
        assert_allclose(reg.outer_height.to_value("arcsec"), 80)
        assert_allclose(reg.angle.to_value("deg"), 0)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, EllipseAnnulusPixelRegion)


class TestRectangleAnnulusPixelRegion(BaseTestPixelRegion):
    reg = RectangleAnnulusPixelRegion(PixCoord(3, 4), 2, 5, 5, 8)
    sample_box = [1, 6, 0, 9]
    inside = [(3, 7)]
    outside = [(3, 4)]
    expected_area = 30
    expected_repr = ('<RectangleAnnulusPixelRegion(PixCoord(x=3, y=4), '
                     'inner width=2, inner height=5, outer width=5, '
                     'outer height=8, angle=0.0 deg)>')
    expected_str = ('Region: RectangleAnnulusPixelRegion\ncenter: PixCoord(x=3, y=4)'
                    '\ninner width: 2\ninner height: 5\nouter width: 5\n'
                    'outer height: 8\nangle: 0.0 deg')

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_width, 2)
        assert_quantity_allclose(self.reg.inner_height, 5)
        assert_quantity_allclose(self.reg.outer_width, 5)
        assert_quantity_allclose(self.reg.outer_height, 8)

    def test_copy(self):
        reg = self.reg.copy()
        assert_quantity_allclose(reg.center.x, 3)
        assert_quantity_allclose(reg.center.y, 4)
        assert_quantity_allclose(reg.inner_width, 2)
        assert_quantity_allclose(reg.inner_height, 5)
        assert_quantity_allclose(reg.outer_width, 5)
        assert_quantity_allclose(reg.outer_height, 8)

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, RectangleAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value("deg"), 90)


class TestRectangleAnnulusSkyRegion(BaseTestSkyRegion):
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    reg = RectangleAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                    50 * u.arcsec, 80 * u.arcsec)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<RectangleAnnulusSkyRegion(<SkyCoord (ICRS): (ra, dec) in '
                     'deg\n    (3., 4.)>, inner width=20.0 arcsec, '
                     'inner height=50.0 arcsec, outer width=50.0 arcsec, '
                     'outer height=80.0 arcsec, angle=0.0 deg)>')
    expected_str = ('Region: RectangleAnnulusSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>\ninner width: 20.0 arcsec\n'
                    'inner height: 50.0 arcsec\nouter width: 50.0 arcsec\n'
                    'outer height: 80.0 arcsec\nangle: 0.0 deg')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_width, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.inner_height, 50 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_width.to_value("arcsec"), 20)
        assert_allclose(reg.inner_height.to_value("arcsec"), 50)
        assert_allclose(reg.outer_width.to_value("arcsec"), 50)
        assert_allclose(reg.outer_height.to_value("arcsec"), 80)
        assert_allclose(reg.angle.to_value("deg"), 0)
        assert reg.visual == {}
        assert reg.meta == {}

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, RectangleAnnulusPixelRegion)
