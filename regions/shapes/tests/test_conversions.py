# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for to_pixel and to_sky coordinate conversions.

Automatic dispatch uses the Jacobian method for distorted WCS
(e.g., SIP) and the offset-based method for non-distorted WCS.
"""

import astropy.units as u
import pytest
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_GWCS
from regions.core import PixCoord
from regions.shapes.annulus import (CircleAnnulusPixelRegion,
                                    CircleAnnulusSkyRegion,
                                    EllipseAnnulusPixelRegion,
                                    EllipseAnnulusSkyRegion,
                                    RectangleAnnulusPixelRegion,
                                    RectangleAnnulusSkyRegion)
from regions.shapes.circle import CirclePixelRegion, CircleSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion
from regions.tests.helpers import WCS_CENTER as CENTER

# -------------------------------------------------------------------
# Region factories: (sky_region, pixel_region) pairs for each shape
# -------------------------------------------------------------------
# Each entry: (sky_cls, sky_kwargs, pix_cls, pix_kwargs,
#              size_attrs, has_angle)
#
# size_attrs: attribute names compared in roundtrip assertions
# has_angle:  whether the region has a rotation angle
DIRECTED_REGIONS = [
    pytest.param(
        EllipseSkyRegion,
        dict(center=CENTER, width=200 * u.arcsec,
             height=100 * u.arcsec, angle=30 * u.deg),
        EllipsePixelRegion,
        dict(center=PixCoord(10, 10), width=5, height=3,
             angle=45 * u.deg),
        ('width', 'height'),
        True,
        id='ellipse',
    ),
    pytest.param(
        RectangleSkyRegion,
        dict(center=CENTER, width=200 * u.arcsec,
             height=100 * u.arcsec, angle=30 * u.deg),
        RectanglePixelRegion,
        dict(center=PixCoord(10, 10), width=5, height=3,
             angle=45 * u.deg),
        ('width', 'height'),
        True,
        id='rectangle',
    ),
    pytest.param(
        EllipseAnnulusSkyRegion,
        dict(center=CENTER,
             inner_width=100 * u.arcsec, outer_width=200 * u.arcsec,
             inner_height=50 * u.arcsec, outer_height=100 * u.arcsec,
             angle=30 * u.deg),
        EllipseAnnulusPixelRegion,
        dict(center=PixCoord(10, 10),
             inner_width=3, outer_width=5,
             inner_height=2, outer_height=4,
             angle=45 * u.deg),
        ('inner_width', 'outer_width', 'inner_height', 'outer_height'),
        True,
        id='ellipse_annulus',
    ),
    pytest.param(
        RectangleAnnulusSkyRegion,
        dict(center=CENTER,
             inner_width=100 * u.arcsec, outer_width=200 * u.arcsec,
             inner_height=50 * u.arcsec, outer_height=100 * u.arcsec,
             angle=30 * u.deg),
        RectangleAnnulusPixelRegion,
        dict(center=PixCoord(10, 10),
             inner_width=3, outer_width=5,
             inner_height=2, outer_height=4,
             angle=45 * u.deg),
        ('inner_width', 'outer_width', 'inner_height', 'outer_height'),
        True,
        id='rectangle_annulus',
    ),
]

CIRCULAR_REGIONS = [
    pytest.param(
        CircleSkyRegion,
        dict(center=CENTER, radius=100 * u.arcsec),
        CirclePixelRegion,
        dict(center=PixCoord(10, 10), radius=3),
        ('radius',),
        False,
        id='circle',
    ),
    pytest.param(
        CircleAnnulusSkyRegion,
        dict(center=CENTER,
             inner_radius=50 * u.arcsec, outer_radius=100 * u.arcsec),
        CircleAnnulusPixelRegion,
        dict(center=PixCoord(10, 10),
             inner_radius=2, outer_radius=4),
        ('inner_radius', 'outer_radius'),
        False,
        id='circle_annulus',
    ),
]

ALL_REGIONS = DIRECTED_REGIONS + CIRCULAR_REGIONS


class TestSkyToPixel:
    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        ALL_REGIONS)
    def test_simple_wcs(self, simple_wcs, sky_cls, sky_kw, pix_cls,
                        pix_kw, size_attrs, has_angle):
        pix_reg = sky_cls(**sky_kw).to_pixel(simple_wcs)
        assert isinstance(pix_reg, pix_cls)
        for attr in size_attrs:
            assert getattr(pix_reg, attr) > 0

    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        DIRECTED_REGIONS)
    def test_rotated_wcs(self, rotated_wcs, sky_cls, sky_kw, pix_cls,
                         pix_kw, size_attrs, has_angle):
        pix_reg = sky_cls(**sky_kw).to_pixel(rotated_wcs)
        assert isinstance(pix_reg, pix_cls)
        for attr in size_attrs:
            assert getattr(pix_reg, attr) > 0

    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        ALL_REGIONS)
    def test_sip_wcs(self, sip_wcs, sky_cls, sky_kw, pix_cls,
                     pix_kw, size_attrs, has_angle):
        """
        Test that the Jacobian path is automatically used for distorted
        WCS.
        """
        pix_reg = sky_cls(**sky_kw).to_pixel(sip_wcs)
        assert isinstance(pix_reg, pix_cls)
        for attr in size_attrs:
            assert getattr(pix_reg, attr) > 0


class TestPixelToSky:
    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        ALL_REGIONS)
    def test_simple_wcs(self, simple_wcs, sky_cls, sky_kw, pix_cls,
                        pix_kw, size_attrs, has_angle):
        sky_reg = pix_cls(**pix_kw).to_sky(simple_wcs)
        assert isinstance(sky_reg, sky_cls)
        for attr in size_attrs:
            assert getattr(sky_reg, attr) > 0 * u.arcsec


class TestRoundtripSkyPixelSky:
    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        ALL_REGIONS)
    def test_simple_wcs(self, simple_wcs, sky_cls, sky_kw, pix_cls,
                        pix_kw, size_attrs, has_angle):
        sky_reg = sky_cls(**sky_kw)
        sky_rt = sky_reg.to_pixel(simple_wcs).to_sky(simple_wcs)
        assert isinstance(sky_rt, sky_cls)
        for attr in size_attrs:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky_reg, attr))
        if has_angle:
            assert_quantity_allclose(sky_rt.angle, sky_reg.angle)

    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        ALL_REGIONS)
    def test_sip_wcs(self, sip_wcs, sky_cls, sky_kw, pix_cls,
                     pix_kw, size_attrs, has_angle):
        """
        Test the roundtrip with distorted WCS (tests Jacobian path).
        """
        sky_reg = sky_cls(**sky_kw)
        sky_rt = sky_reg.to_pixel(sip_wcs).to_sky(sip_wcs)
        for attr in size_attrs:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky_reg, attr))
        if has_angle:
            assert_quantity_allclose(sky_rt.angle, sky_reg.angle)


class TestRoundtripPixelSkyPixel:
    @pytest.mark.parametrize(
        'sky_cls, sky_kw, pix_cls, pix_kw, size_attrs, has_angle',
        DIRECTED_REGIONS)
    def test_simple_wcs(self, simple_wcs, sky_cls, sky_kw, pix_cls,
                        pix_kw, size_attrs, has_angle):
        pix_reg = pix_cls(**pix_kw)
        pix_rt = pix_reg.to_sky(simple_wcs).to_pixel(simple_wcs)
        assert_allclose(pix_rt.center.x, pix_reg.center.x)
        assert_allclose(pix_rt.center.y, pix_reg.center.y)
        for attr in size_attrs:
            assert_allclose(getattr(pix_rt, attr), getattr(pix_reg, attr))
        if has_angle:
            assert_quantity_allclose(pix_rt.angle, pix_reg.angle)


# Build GWCS test region params with smaller sky sizes (0.05 deg/pix)
GWCS_SKY_REGIONS = [
    pytest.param(
        CircleSkyRegion,
        dict(center=CENTER, radius=5 * u.arcsec),
        ('radius',), False,
        id='circle',
    ),
    pytest.param(
        EllipseSkyRegion,
        dict(center=CENTER, width=10 * u.arcsec,
             height=5 * u.arcsec, angle=30 * u.deg),
        ('width', 'height'), True,
        id='ellipse',
    ),
    pytest.param(
        RectangleSkyRegion,
        dict(center=CENTER, width=10 * u.arcsec,
             height=5 * u.arcsec, angle=30 * u.deg),
        ('width', 'height'), True,
        id='rectangle',
    ),
    pytest.param(
        CircleAnnulusSkyRegion,
        dict(center=CENTER,
             inner_radius=3 * u.arcsec, outer_radius=6 * u.arcsec),
        ('inner_radius', 'outer_radius'), False,
        id='circle_annulus',
    ),
    pytest.param(
        EllipseAnnulusSkyRegion,
        dict(center=CENTER,
             inner_width=6 * u.arcsec, outer_width=12 * u.arcsec,
             inner_height=3 * u.arcsec, outer_height=6 * u.arcsec,
             angle=30 * u.deg),
        ('inner_width', 'outer_width', 'inner_height', 'outer_height'),
        True,
        id='ellipse_annulus',
    ),
    pytest.param(
        RectangleAnnulusSkyRegion,
        dict(center=CENTER,
             inner_width=6 * u.arcsec, outer_width=12 * u.arcsec,
             inner_height=3 * u.arcsec, outer_height=6 * u.arcsec,
             angle=30 * u.deg),
        ('inner_width', 'outer_width', 'inner_height', 'outer_height'),
        True,
        id='rectangle_annulus',
    ),
]


@pytest.mark.skipif(not HAS_GWCS, reason='gwcs is required')
class TestGWCSConversion:
    """
    Tests for coordinate conversions with GWCS objects.

    GWCS objects are always treated as distorted (Jacobian path) because
    they do not have the ``has_distortion`` attribute used by astropy
    WCS.
    """

    def test_gwcs_has_distortion(self, gwcs_obj):
        assert getattr(gwcs_obj, 'has_distortion', True)

    @pytest.mark.parametrize(
        'sky_cls, sky_kw, size_attrs, has_angle', GWCS_SKY_REGIONS)
    def test_roundtrip(self, gwcs_obj, sky_cls, sky_kw, size_attrs,
                       has_angle):
        sky_reg = sky_cls(**sky_kw)
        sky_rt = sky_reg.to_pixel(gwcs_obj).to_sky(gwcs_obj)
        assert isinstance(sky_rt, sky_cls)
        for attr in size_attrs:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky_reg, attr))
        if has_angle:
            assert_quantity_allclose(sky_rt.angle, sky_reg.angle)


class TestDistortionDetection:
    """
    Test WCS distortion detection used for automatic method dispatch.
    """

    def test_simple_wcs_has_no_distortion(self, simple_wcs):
        assert not getattr(simple_wcs, 'has_distortion', True)

    def test_sip_wcs_has_distortion(self, sip_wcs):
        assert getattr(sip_wcs, 'has_distortion', True)

    @pytest.mark.skipif(not HAS_GWCS, reason='gwcs is required')
    def test_gwcs_has_distortion(self, gwcs_obj):
        assert getattr(gwcs_obj, 'has_distortion', True)
