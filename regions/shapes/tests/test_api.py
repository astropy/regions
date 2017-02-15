# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The tests in this file simply check what functionality is currently
implemented and doesn't check anything about correctness.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import itertools
from astropy.tests.helper import pytest
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.wcs import WCS
from ...core.mask import Mask
from ...core.core import Region, SkyRegion, PixelRegion
from ...core.pixcoord import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from .utils import HAS_SHAPELY

PIXEL_REGIONS = [
    CirclePixelRegion(PixCoord(3, 4), radius=5),
    EllipsePixelRegion(PixCoord(3, 4), major=7, minor=5, angle=3 * u.deg),
    PolygonPixelRegion(PixCoord([1, 4, 3], [2, 4, 4])),
    RectanglePixelRegion(PixCoord(6, 5), width=3, height=5)
]

SKY_REGIONS = [
    CircleSkyRegion(ICRS(3 * u.deg, 4 * u.deg), radius=5 * u.deg),
    EllipseSkyRegion(ICRS(3 * u.deg, 4 * u.deg), major=7 * u.deg, minor=5 * u.deg, angle=3 * u.deg),
    PolygonSkyRegion(ICRS([1, 4, 3] * u.deg, [2, 4, 4] * u.deg)),
    RectangleSkyRegion(ICRS(6 * u.deg, 5 * u.deg), width=3 * u.deg, height=5 * u.deg)
]

SKYPIX_MODES = ['local', 'affine', 'full']
MASK_MODES = ['center', 'exact', 'subpixels']
COMMON_WCS = WCS(naxis=2)
COMMON_WCS.wcs.ctype = 'RA---TAN', 'DEC--TAN'


def ids_func(arg):
    if isinstance(arg, Region):
        return arg.__class__.__name__
    else:
        return str(arg)


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_in(region):
    try:
        PixCoord(1, 1) in region
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_area(region):
    try:
        area = region.area
        assert not isinstance(area, u.Quantity)
    except AttributeError:
        pytest.xfail()


@pytest.mark.parametrize(('region', 'mode'), itertools.product(PIXEL_REGIONS, SKYPIX_MODES), ids=ids_func)
def test_pix_to_sky(region, mode):
    try:
        sky_region = region.to_sky(COMMON_WCS, mode=mode)
        assert isinstance(sky_region, SkyRegion)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.skipif('not HAS_SHAPELY')
@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_to_shapely(region):
    try:
        from shapely.geometry.base import BaseGeometry
        shape = region.to_shapely()
        assert isinstance(shape, BaseGeometry)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize(('region', 'mode'), itertools.product(PIXEL_REGIONS, MASK_MODES), ids=ids_func)
def test_pix_to_mask(region, mode):
    try:
        mask = region.to_mask(mode=mode)
        assert isinstance(mask, Mask)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in(region):
    try:
        ICRS(1 * u.deg, 1 * u.deg) in region
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_area(region):
    try:
        area = region.area
        assert isinstance(area, u.Quantity)
    except AttributeError:
        pytest.xfail()


@pytest.mark.parametrize(
    ('region', 'mode'),
    itertools.product(SKY_REGIONS, SKYPIX_MODES),
    ids=ids_func,
)
def test_sky_to_pix(region, mode):
    try:
        pix_region = region.to_pixel(mode=mode, wcs=COMMON_WCS)
        assert isinstance(pix_region, PixelRegion)
    except NotImplementedError:
        pytest.xfail()
