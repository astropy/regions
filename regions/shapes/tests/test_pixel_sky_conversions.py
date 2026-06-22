# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for region pixel-to-sky and sky-to-pixel coordinate conversions.

This module covers:

* roundtrip conversions for all pixel and sky region classes through
  simple, rotated, distorted (SIP), sheared, flipped-parity, and
  generalized (`gwcs.wcs.WCS`) WCS objects, exercising the SVD/Jacobian
  shape conversion (the SIP and gWCS cases in particular force the
  distortion (Jacobian) code path of the internal ``_wcs_helpers``)

* the rotation-angle convention for directed (elliptical and
  rectangular) regions

Automatic dispatch uses the Jacobian method for distorted WCS (e.g.,
SIP) and the offset-based method for non-distorted WCS.

The regions rotation-angle convention measures the angle of the
region's width axis from the longitude (RA) axis (sky regions) or the
positive ``x`` axis (pixel regions), counterclockwise. This differs
from the photutils convention (position angle from North), so for an
axis-aligned, isotropic WCS the regions sky angle and pixel angle are
numerically *equal* after conversion (no 90 deg offset).
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.io.fits import Header
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_GWCS
from regions.core import PixCoord
from regions.core.metadata import RegionVisual
from regions.shapes.annulus import (CircleAnnulusPixelRegion,
                                    CircleAnnulusSkyRegion,
                                    EllipseAnnulusPixelRegion,
                                    EllipseAnnulusSkyRegion,
                                    RectangleAnnulusPixelRegion,
                                    RectangleAnnulusSkyRegion)
from regions.shapes.circle import CirclePixelRegion, CircleSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion
from regions.shapes.text import TextPixelRegion, TextSkyRegion
from regions.tests.helpers import WCS_CENTER as CENTER
from regions.tests.helpers import make_gwcs

# Module constants
PIX_CENTER = (10, 10)

# Rotation angles (deg) used to pin the sky/pixel angle convention.
ANGLES_DEG = [0, 30, 60, 135, 200, 315]

# Rotation angles (deg) used for the sheared-WCS and degenerate-shape
# roundtrip tests.
SHEAR_ANGLES_DEG = [0.0, 30.0, 75.0]

# Each case bundles a matching (sky, pixel) region class pair with their
# construction kwargs. ``size_attrs`` lists the shape attributes
# compared in roundtrip assertions and ``has_angle`` flags whether the
# region carries an ``angle`` rotation.
_REGION_CASES = [
    {'id': 'ellipse',
     'sky_cls': EllipseSkyRegion,
     'sky_kw': {'width': 200 * u.arcsec, 'height': 100 * u.arcsec,
                'angle': 30 * u.deg},
     'pix_cls': EllipsePixelRegion,
     'pix_kw': {'width': 5, 'height': 3, 'angle': 45 * u.deg},
     'size_attrs': ('width', 'height'),
     'has_angle': True},
    {'id': 'rectangle',
     'sky_cls': RectangleSkyRegion,
     'sky_kw': {'width': 200 * u.arcsec, 'height': 100 * u.arcsec,
                'angle': 30 * u.deg},
     'pix_cls': RectanglePixelRegion,
     'pix_kw': {'width': 5, 'height': 3, 'angle': 45 * u.deg},
     'size_attrs': ('width', 'height'),
     'has_angle': True},
    {'id': 'ellipse_annulus',
     'sky_cls': EllipseAnnulusSkyRegion,
     'sky_kw': {'inner_width': 100 * u.arcsec, 'outer_width': 200 * u.arcsec,
                'inner_height': 50 * u.arcsec, 'outer_height': 100 * u.arcsec,
                'angle': 30 * u.deg},
     'pix_cls': EllipseAnnulusPixelRegion,
     'pix_kw': {'inner_width': 3, 'outer_width': 6,
                'inner_height': 2, 'outer_height': 4, 'angle': 45 * u.deg},
     'size_attrs': ('inner_width', 'outer_width', 'inner_height',
                    'outer_height'),
     'has_angle': True},
    {'id': 'rectangle_annulus',
     'sky_cls': RectangleAnnulusSkyRegion,
     'sky_kw': {'inner_width': 100 * u.arcsec, 'outer_width': 200 * u.arcsec,
                'inner_height': 50 * u.arcsec, 'outer_height': 100 * u.arcsec,
                'angle': 30 * u.deg},
     'pix_cls': RectangleAnnulusPixelRegion,
     'pix_kw': {'inner_width': 3, 'outer_width': 6,
                'inner_height': 2, 'outer_height': 4, 'angle': 45 * u.deg},
     'size_attrs': ('inner_width', 'outer_width', 'inner_height',
                    'outer_height'),
     'has_angle': True},
    {'id': 'circle',
     'sky_cls': CircleSkyRegion,
     'sky_kw': {'radius': 100 * u.arcsec},
     'pix_cls': CirclePixelRegion,
     'pix_kw': {'radius': 3},
     'size_attrs': ('radius',),
     'has_angle': False},
    {'id': 'circle_annulus',
     'sky_cls': CircleAnnulusSkyRegion,
     'sky_kw': {'inner_radius': 50 * u.arcsec, 'outer_radius': 100 * u.arcsec},
     'pix_cls': CircleAnnulusPixelRegion,
     'pix_kw': {'inner_radius': 2, 'outer_radius': 4},
     'size_attrs': ('inner_radius', 'outer_radius'),
     'has_angle': False},
]

# All region cases and the directed (rotatable) subset, wrapped as
# pytest params with readable ids.
REGION_CASES = [pytest.param(case, id=case['id'])
                for case in _REGION_CASES]
DIRECTED_CASES = [pytest.param(case, id=case['id'])
                  for case in _REGION_CASES if case['has_angle']]


def _angle_diff_deg(actual, desired):
    """
    Signed angular difference in degrees, wrapped to (-180, 180].
    """
    diff = (actual.to_value(u.deg) - desired.to_value(u.deg) + 180) % 360
    return diff - 180


def _build_pair(case, angle):
    """
    Build a matching (sky_region, pixel_region) pair from a case,
    overriding the rotation ``angle`` of both regions.
    """
    sky = case['sky_cls'](CENTER, **{**case['sky_kw'], 'angle': angle})
    pix = case['pix_cls'](PixCoord(*PIX_CENTER),
                          **{**case['pix_kw'], 'angle': angle})
    return sky, pix


def _make_sip_wcs(ra_deg=100.0, dec_deg=30.0):
    """
    Build a small TAN-SIP WCS centered at the given (RA, Dec).

    The SIP terms are tiny but nonzero, which forces the Jacobian
    (distortion) code path to be exercised by ``to_pixel``/``to_sky``.
    """
    header = Header()
    header['NAXIS'] = 2
    header['NAXIS1'] = 21
    header['NAXIS2'] = 21
    header['CRPIX1'] = 10.5
    header['CRPIX2'] = 10.5
    header['CRVAL1'] = ra_deg
    header['CRVAL2'] = dec_deg
    header['CTYPE1'] = 'RA---TAN-SIP'
    header['CTYPE2'] = 'DEC--TAN-SIP'
    cdelt = 0.1 / 3600.0
    header['CD1_1'] = -cdelt
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = cdelt
    header['A_ORDER'] = 2
    header['A_2_0'] = 1e-7
    header['A_0_2'] = 1e-7
    header['B_ORDER'] = 2
    header['B_2_0'] = 1e-7
    header['B_0_2'] = 1e-7

    return WCS(header)


@pytest.fixture
def simple_wcs():
    """
    A simple axis-aligned TAN WCS (RA increases to the left, equal pixel
    scales, North = +y).
    """
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [10.5, 10.5]
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.cdelt = [-0.1 / 3600, 0.1 / 3600]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcs


@pytest.fixture
def rotated_wcs():
    """
    A 25 deg rotated TAN WCS.
    """
    rotation_deg = 25.0
    cdelt = 0.1 / 3600
    rad = np.radians(rotation_deg)
    cos_a = np.cos(rad)
    sin_a = np.sin(rad)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [10.5, 10.5]
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.cd = [[-cdelt * cos_a, cdelt * sin_a],
                  [cdelt * sin_a, cdelt * cos_a]]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcs


@pytest.fixture
def sheared_wcs():
    """
    A sheared TAN WCS where the pixel x and y axes are 80 deg apart on
    the sky (10 deg of shear from a perpendicular axis-aligned grid).
    """
    cdelt = 0.1 / 3600
    a = np.radians(80.0)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [10.5, 10.5]
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.cd = [[-cdelt, cdelt * np.cos(a)],
                  [0.0, cdelt * np.sin(a)]]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcs


@pytest.fixture
def sip_wcs():
    """
    A TAN-SIP WCS with small distortion terms centered at CENTER.
    """
    return _make_sip_wcs()


@pytest.fixture
def flipped_wcs():
    """
    A flipped-parity TAN WCS (North down, East left).

    Both CDELT values are negative, so the pixel scale matrix has
    a positive determinant (parity = +1), opposite the standard
    astronomical convention. North (increasing Dec) points along -y and
    East (increasing RA) points along -x.
    """
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [10.5, 10.5]
    wcs.wcs.crval = [CENTER.ra.deg, CENTER.dec.deg]
    wcs.wcs.cdelt = [-0.1 / 3600, -0.1 / 3600]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcs


class TestSkyToPixel:
    """
    Converting a sky region to a pixel region must return the matching
    pixel class with positive shape parameters.
    """

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_simple_wcs(self, simple_wcs, case):
        pix = case['sky_cls'](CENTER, **case['sky_kw']).to_pixel(simple_wcs)
        assert isinstance(pix, case['pix_cls'])
        for attr in case['size_attrs']:
            assert getattr(pix, attr) > 0

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_rotated_wcs(self, rotated_wcs, case):
        pix = case['sky_cls'](CENTER, **case['sky_kw']).to_pixel(rotated_wcs)
        assert isinstance(pix, case['pix_cls'])
        for attr in case['size_attrs']:
            assert getattr(pix, attr) > 0

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_sip_wcs(self, sip_wcs, case):
        """
        The Jacobian (distortion) path is used automatically for a SIP
        WCS.
        """
        pix = case['sky_cls'](CENTER, **case['sky_kw']).to_pixel(sip_wcs)
        assert isinstance(pix, case['pix_cls'])
        for attr in case['size_attrs']:
            assert getattr(pix, attr) > 0


class TestPixelToSky:
    """
    Converting a pixel region to a sky region must return the matching
    sky class with positive angular shape parameters.
    """

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_simple_wcs(self, simple_wcs, case):
        sky = case['pix_cls'](PixCoord(*PIX_CENTER),
                              **case['pix_kw']).to_sky(simple_wcs)
        assert isinstance(sky, case['sky_cls'])
        for attr in case['size_attrs']:
            assert getattr(sky, attr) > 0 * u.arcsec

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_sip_wcs(self, sip_wcs, case):
        sky = case['pix_cls'](PixCoord(*PIX_CENTER),
                              **case['pix_kw']).to_sky(sip_wcs)
        assert isinstance(sky, case['sky_cls'])
        for attr in case['size_attrs']:
            assert getattr(sky, attr) > 0 * u.arcsec


class TestRoundtripSkyPixelSky:
    """
    A sky -> pixel -> sky roundtrip must recover the original shape
    parameters and rotation angle.
    """

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_simple_wcs(self, simple_wcs, case):
        sky = case['sky_cls'](CENTER, **case['sky_kw'])
        sky_rt = sky.to_pixel(simple_wcs).to_sky(simple_wcs)
        assert isinstance(sky_rt, case['sky_cls'])
        for attr in case['size_attrs']:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky, attr))
        if case['has_angle']:
            assert_quantity_allclose(sky_rt.angle, sky.angle,
                                     atol=1e-9 * u.deg)

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_rotated_wcs(self, rotated_wcs, case):
        sky = case['sky_cls'](CENTER, **case['sky_kw'])
        sky_rt = sky.to_pixel(rotated_wcs).to_sky(rotated_wcs)
        for attr in case['size_attrs']:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky, attr))
        if case['has_angle']:
            assert_quantity_allclose(sky_rt.angle, sky.angle,
                                     atol=1e-9 * u.deg)

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_sip_wcs(self, sip_wcs, case):
        """
        Roundtrip through a distorted (SIP) WCS, exercising the Jacobian
        code path.
        """
        sky = case['sky_cls'](CENTER, **case['sky_kw'])
        sky_rt = sky.to_pixel(sip_wcs).to_sky(sip_wcs)
        for attr in case['size_attrs']:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky, attr), rtol=1e-3)
        if case['has_angle']:
            assert_quantity_allclose(sky_rt.angle, sky.angle,
                                     atol=1e-3 * u.deg)


class TestRoundtripPixelSkyPixel:
    """
    A pixel -> sky -> pixel roundtrip must recover the original center,
    shape parameters, and rotation angle.
    """

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_simple_wcs(self, simple_wcs, case):
        pix = case['pix_cls'](PixCoord(*PIX_CENTER), **case['pix_kw'])
        pix_rt = pix.to_sky(simple_wcs).to_pixel(simple_wcs)
        assert_allclose(pix_rt.center.x, pix.center.x)
        assert_allclose(pix_rt.center.y, pix.center.y)
        for attr in case['size_attrs']:
            assert_allclose(getattr(pix_rt, attr), getattr(pix, attr))
        if case['has_angle']:
            assert_quantity_allclose(pix_rt.angle, pix.angle,
                                     atol=1e-9 * u.deg)

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_rotated_wcs(self, rotated_wcs, case):
        pix = case['pix_cls'](PixCoord(*PIX_CENTER), **case['pix_kw'])
        pix_rt = pix.to_sky(rotated_wcs).to_pixel(rotated_wcs)
        assert_allclose(pix_rt.center.x, pix.center.x)
        assert_allclose(pix_rt.center.y, pix.center.y)
        for attr in case['size_attrs']:
            assert_allclose(getattr(pix_rt, attr), getattr(pix, attr))
        if case['has_angle']:
            assert_quantity_allclose(pix_rt.angle, pix.angle,
                                     atol=1e-9 * u.deg)

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_sip_wcs(self, sip_wcs, case):
        pix = case['pix_cls'](PixCoord(*PIX_CENTER), **case['pix_kw'])
        pix_rt = pix.to_sky(sip_wcs).to_pixel(sip_wcs)
        assert_allclose(pix_rt.center.x, pix.center.x, atol=1e-6)
        assert_allclose(pix_rt.center.y, pix.center.y, atol=1e-6)
        for attr in case['size_attrs']:
            assert_allclose(getattr(pix_rt, attr), getattr(pix, attr),
                            rtol=1e-3)
        if case['has_angle']:
            assert_quantity_allclose(pix_rt.angle, pix.angle,
                                     atol=1e-3 * u.deg)


@pytest.mark.skipif(not HAS_GWCS, reason='gwcs is required')
class TestGWCSRoundtrip:
    """
    Roundtrip conversions through a generalized WCS (`gwcs.wcs.WCS`).

    A gWCS has no ``has_distortion`` attribute, so the Jacobian path is
    always used.
    """

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_sky_pixel_sky(self, case):
        wcs = make_gwcs()
        sky = case['sky_cls'](CENTER, **case['sky_kw'])
        sky_rt = sky.to_pixel(wcs).to_sky(wcs)
        assert isinstance(sky_rt, case['sky_cls'])
        for attr in case['size_attrs']:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky, attr), rtol=1e-5)
        if case['has_angle']:
            assert_quantity_allclose(sky_rt.angle, sky.angle,
                                     atol=1e-4 * u.deg)

    @pytest.mark.parametrize('case', REGION_CASES)
    def test_pixel_sky_pixel(self, case):
        wcs = make_gwcs()
        pix = case['pix_cls'](PixCoord(*PIX_CENTER), **case['pix_kw'])
        pix_rt = pix.to_sky(wcs).to_pixel(wcs)
        assert_allclose(pix_rt.center.x, pix.center.x, rtol=1e-5)
        assert_allclose(pix_rt.center.y, pix.center.y, rtol=1e-5)
        for attr in case['size_attrs']:
            assert_allclose(getattr(pix_rt, attr), getattr(pix, attr),
                            rtol=1e-5)
        if case['has_angle']:
            assert_quantity_allclose(pix_rt.angle, pix.angle,
                                     atol=1e-4 * u.deg)


class TestSkyPixelAngleConvention:
    """
    Verify the sky/pixel angle convention on an axis-aligned WCS.

    The regions convention measures the rotation angle of the region's
    width axis from the longitude (RA) axis (sky regions) or the
    positive ``x`` axis (pixel regions), counterclockwise. Internally
    the WCS helpers use the position-angle (from North) convention; the
    regions shape methods apply a 90 deg offset to bridge the two.

    For a standard axis-aligned WCS (RA increasing to the left, equal
    isotropic pixel scale, North = +y), the regions sky angle and pixel
    angle must be numerically equal after conversion. These tests pin
    that convention down so that any future refactor of the internal
    ``_wcs_helpers`` cannot silently re-introduce a symmetric 90 deg
    offset bug, which is invisible to roundtrip-only tests because the
    offset cancels in sky -> pixel -> sky (and vice versa).
    """

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', ANGLES_DEG)
    def test_sky_to_pixel_angle(self, simple_wcs, case, angle_deg):
        """
        For an axis-aligned WCS, ``sky.to_pixel`` produces a pixel region
        with the same numerical rotation angle as the sky region.
        """
        sky, _ = _build_pair(case, angle_deg * u.deg)
        pix = sky.to_pixel(simple_wcs)
        assert abs(_angle_diff_deg(pix.angle, angle_deg * u.deg)) < 1e-3

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', ANGLES_DEG)
    def test_pixel_to_sky_angle(self, simple_wcs, case, angle_deg):
        """
        For an axis-aligned WCS, ``pixel.to_sky`` produces a sky region
        with the same numerical rotation angle as the pixel region.
        """
        _, pix = _build_pair(case, angle_deg * u.deg)
        sky = pix.to_sky(simple_wcs)
        assert abs(_angle_diff_deg(sky.angle, angle_deg * u.deg)) < 1e-3

    def test_text_region_rotation_sky_to_pixel(self, simple_wcs):
        """
        ``TextSkyRegion.to_pixel`` must preserve the visual rotation
        angle (in degrees) for an axis-aligned WCS.
        """
        sky = TextSkyRegion(CENTER, 'foo', visual=RegionVisual(rotation=30.0))
        pix = sky.to_pixel(simple_wcs)
        diff = ((pix.visual['rotation'] - 30.0 + 180) % 360) - 180
        assert abs(diff) < 1e-3

    def test_text_region_rotation_pixel_to_sky(self, simple_wcs):
        """
        ``TextPixelRegion.to_sky`` must preserve the visual rotation
        angle (in degrees) for an axis-aligned WCS.
        """
        pix = TextPixelRegion(PixCoord(9.5, 9.5), 'foo',
                              visual=RegionVisual(rotation=30.0))
        sky = pix.to_sky(simple_wcs)
        diff = ((sky.visual['rotation'] - 30.0 + 180) % 360) - 180
        assert abs(diff) < 1e-3


class TestShearedWCSRoundtrip:
    """
    Verify that directed regions round-trip exactly through a sheared
    WCS, where the pixel x and y axes are not perpendicular on the sky.
    """

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    def test_sky_pixel_sky(self, sheared_wcs, case, angle_deg):
        sky, _ = _build_pair(case, angle_deg * u.deg)
        sky_rt = sky.to_pixel(sheared_wcs).to_sky(sheared_wcs)
        assert abs(_angle_diff_deg(sky_rt.angle, angle_deg * u.deg)) < 1e-6
        for attr in case['size_attrs']:
            assert_quantity_allclose(getattr(sky_rt, attr),
                                     getattr(sky, attr), rtol=1e-6)

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    def test_pixel_sky_pixel(self, sheared_wcs, case, angle_deg):
        _, pix = _build_pair(case, angle_deg * u.deg)
        pix_rt = pix.to_sky(sheared_wcs).to_pixel(sheared_wcs)
        assert abs(_angle_diff_deg(pix_rt.angle, angle_deg * u.deg)) < 1e-6
        for attr in case['size_attrs']:
            assert_allclose(getattr(pix_rt, attr), getattr(pix, attr),
                            rtol=1e-6)

    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    def test_text_sky_pixel_sky_rotation(self, sheared_wcs, angle_deg):
        """
        Text region rotation must round-trip through a sheared WCS. The
        SVD path preserves the input theta direction when the input
        ellipse is circular (no width/height), via the degenerate
        singular-value fallback in the shape helper.
        """
        sky = TextSkyRegion(CENTER, 'foo',
                            visual=RegionVisual(rotation=angle_deg))
        rt = sky.to_pixel(sheared_wcs).to_sky(sheared_wcs)
        diff = ((rt.visual['rotation'] - angle_deg + 180) % 360) - 180
        assert abs(diff) < 1e-6

    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    def test_text_pixel_sky_pixel_rotation(self, sheared_wcs, angle_deg):
        pix = TextPixelRegion(PixCoord(*PIX_CENTER), 'foo',
                              visual=RegionVisual(rotation=angle_deg))
        rt = pix.to_sky(sheared_wcs).to_pixel(sheared_wcs)
        diff = ((rt.visual['rotation'] - angle_deg + 180) % 360) - 180
        assert abs(diff) < 1e-6


class TestCircularInputAnglePreserved:
    """
    Verify that the SVD path preserves the input rotation angle when the
    input rectangular/elliptical region is shape-degenerate (width ==
    height): the converted angle must round-trip exactly.

    For a circular shape the SVD principal axis is otherwise arbitrary,
    so the helper falls back to the mapped width semi-axis direction to
    preserve orientation.
    """

    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    @pytest.mark.parametrize(
        'sky_cls', [RectangleSkyRegion, EllipseSkyRegion])
    def test_square_sky_roundtrip(self, rotated_wcs, sky_cls, angle_deg):
        sky = sky_cls(CENTER, width=200 * u.arcsec, height=200 * u.arcsec,
                      angle=angle_deg * u.deg)
        sky_rt = sky.to_pixel(rotated_wcs).to_sky(rotated_wcs)
        assert abs(_angle_diff_deg(sky_rt.angle, angle_deg * u.deg)) < 1e-6

    @pytest.mark.parametrize('angle_deg', SHEAR_ANGLES_DEG)
    @pytest.mark.parametrize(
        'pix_cls', [RectanglePixelRegion, EllipsePixelRegion])
    def test_square_pixel_roundtrip(self, rotated_wcs, pix_cls, angle_deg):
        pix = pix_cls(PixCoord(*PIX_CENTER), width=4, height=4,
                      angle=angle_deg * u.deg)
        pix_rt = pix.to_sky(rotated_wcs).to_pixel(rotated_wcs)
        assert abs(_angle_diff_deg(pix_rt.angle, angle_deg * u.deg)) < 1e-6


class TestFlippedParityWCS:
    """
    Regression tests for a flipped-parity WCS (North down, East left;
    positive-determinant pixel scale matrix).

    Such a WCS previously produced regions that were mirrored about the
    x-axis. For the flipped axis-aligned WCS, the regions sky and pixel
    rotation angles are related by a mirror, ``pixel.angle == 360 deg -
    sky.angle`` (mod 360), rather than the identity relation that holds
    for a standard-parity axis-aligned WCS. (Regions measures the angle
    from the longitude/x axis, so there is no 90 deg offset as in the
    photutils position-angle convention.)
    """

    def test_flipped_wcs_has_positive_parity(self, flipped_wcs):
        """
        The flipped WCS fixture must have a positive-determinant pixel
        scale matrix (parity = +1).
        """
        assert np.linalg.det(flipped_wcs.pixel_scale_matrix) > 0

    def test_north_is_minus_y(self, flipped_wcs):
        """
        A sky region whose width axis points North (regions sky angle =
        90 deg, i.e. along the Dec axis) maps to a pixel region pointing
        along -y (down), i.e., pixel angle = 270 deg.
        """
        sky = EllipseSkyRegion(CENTER, width=200 * u.arcsec,
                               height=100 * u.arcsec, angle=90 * u.deg)
        pix = sky.to_pixel(flipped_wcs)
        assert abs(_angle_diff_deg(pix.angle, 270 * u.deg)) < 2e-5

    def test_east_is_minus_x(self, flipped_wcs):
        """
        A sky region whose width axis points East (regions sky angle =
        180 deg) maps to a pixel region pointing along -x (left), i.e.,
        pixel angle = 180 deg.
        """
        sky = EllipseSkyRegion(CENTER, width=200 * u.arcsec,
                               height=100 * u.arcsec, angle=180 * u.deg)
        pix = sky.to_pixel(flipped_wcs)
        assert abs(_angle_diff_deg(pix.angle, 180 * u.deg)) < 2e-5

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', ANGLES_DEG)
    def test_sky_to_pixel_angle(self, flipped_wcs, case, angle_deg):
        """
        Verify ``sky.to_pixel`` returns ``360 deg - angle`` (mod 360)
        for the flipped-parity WCS (not the unmirrored value).
        """
        sky, _ = _build_pair(case, angle_deg * u.deg)
        pix = sky.to_pixel(flipped_wcs)
        assert abs(_angle_diff_deg(pix.angle, -angle_deg * u.deg)) < 2e-5

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', ANGLES_DEG)
    def test_pixel_to_sky_angle(self, flipped_wcs, case, angle_deg):
        """
        Verify ``pixel.to_sky`` returns ``360 deg - angle`` (mod 360)
        for the flipped-parity WCS (the inverse of the above relation).
        """
        _, pix = _build_pair(case, angle_deg * u.deg)
        sky = pix.to_sky(flipped_wcs)
        assert abs(_angle_diff_deg(sky.angle, -angle_deg * u.deg)) < 2e-5

    @pytest.mark.parametrize('case', DIRECTED_CASES)
    @pytest.mark.parametrize('angle_deg', ANGLES_DEG)
    def test_sky_pixel_sky_roundtrip_angle(self, flipped_wcs, case, angle_deg):
        """
        Verify ``sky -> pixel -> sky`` preserves the original angle for
        the flipped-parity WCS.
        """
        sky, _ = _build_pair(case, angle_deg * u.deg)
        sky_rt = sky.to_pixel(flipped_wcs).to_sky(flipped_wcs)
        assert abs(_angle_diff_deg(sky_rt.angle, angle_deg * u.deg)) < 2e-5
