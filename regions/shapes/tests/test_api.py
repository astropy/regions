# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The tests in this file simply check what functionality is currently
implemented and doesn't check anything about correctness.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import itertools

from astropy.tests.helper import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from ...core.mask import Mask
from ...core.core import Region, SkyRegion, PixelRegion
from ...core.pixcoord import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from .utils import HAS_SHAPELY  # noqa

PIXEL_REGIONS = [
    CirclePixelRegion(PixCoord(3, 4), radius=5),
    EllipsePixelRegion(PixCoord(3, 4), width=7, height=5, angle=3 * u.deg),
    PolygonPixelRegion(PixCoord([1, 4, 3], [2, 4, 4])),
    RectanglePixelRegion(PixCoord(6, 5), width=3, height=5)]

SKY_REGIONS = [
    CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), radius=5 * u.deg),
    EllipseSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), width=7 * u.deg, height=5 * u.deg, angle=3 * u.deg),
    PolygonSkyRegion(SkyCoord([1, 4, 3] * u.deg, [2, 4, 4] * u.deg)),
    RectangleSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg), width=3 * u.deg, height=5 * u.deg)]

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
    PixCoord(1, 1) in region


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_area(region):
    try:
        area = region.area
        assert not isinstance(area, u.Quantity)
    except ImportError:  # for shapely
        pytest.skip()


@pytest.mark.parametrize(('region'), PIXEL_REGIONS, ids=ids_func)
def test_pix_to_sky(region):
    try:
        sky_region = region.to_sky(COMMON_WCS)
        assert isinstance(sky_region, SkyRegion)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.skipif('not HAS_SHAPELY')
@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_to_shapely(region):
    from shapely.geometry.base import BaseGeometry
    shape = region.to_shapely()
    assert isinstance(shape, BaseGeometry)


@pytest.mark.parametrize(('region', 'mode'),
                         itertools.product(PIXEL_REGIONS, MASK_MODES),
                         ids=ids_func)
def test_pix_to_mask(region, mode):
    try:
        mask = region.to_mask(mode=mode)
        assert isinstance(mask, Mask)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in(region):
    region.contains(SkyCoord(1 * u.deg, 1 * u.deg, frame='icrs'), COMMON_WCS)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_to_pix(region):
    pix_region = region.to_pixel(wcs=COMMON_WCS)
    assert isinstance(pix_region, PixelRegion)
