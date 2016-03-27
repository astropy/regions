# The tests in this file simply check what functionality is currently
# implemented and doesn't check anything about correctness.

import itertools
import pytest
import numpy as np

from astropy import units as u
from astropy.coordinates import ICRS
from astropy.wcs import WCS

from ...core.core import SkyRegion, PixelRegion
from ...core.pixcoord import PixCoord

from ..circle import CirclePixelRegion, CircleSkyRegion
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion

PIXEL_REGIONS = [
    CirclePixelRegion(PixCoord(3, 4), radius=5),
    EllipsePixelRegion(PixCoord(3, 4), minor=5, major=7, angle=3 * u.deg),
    PolygonPixelRegion(PixCoord([1,4,3], [2,4,4])),
    RectanglePixelRegion(PixCoord(6,5), width=3, height=5)
]

SKY_REGIONS = [
    CircleSkyRegion(ICRS(3 * u.deg, 4 * u.deg), radius=5),
    EllipseSkyRegion(ICRS(3 * u.deg, 4 * u.deg), minor=5, major=7, angle=3 * u.deg),
    PolygonSkyRegion(ICRS([1,4,3] * u.deg, [2,4,4] * u.deg)),
    RectangleSkyRegion(ICRS(6 * u.deg,5 * u.deg), width=3, height=5)
]

SKYPIX_MODES = ['local', 'affine', 'full']
MASK_MODES = ['center', 'any', 'all', 'exact']


@pytest.mark.parametrize('region', PIXEL_REGIONS)
def test_pix_in(region):
    PixCoord(1,1) in region


@pytest.mark.parametrize('region', PIXEL_REGIONS)
def test_pix_area(region):
    area = region.area
    assert not isinstance(area, u.Quantity)


@pytest.mark.parametrize(('region', 'mode'), itertools.product(PIXEL_REGIONS, SKYPIX_MODES))
def test_pix_to_sky(region, mode):
    sky_region = region.to_sky(mode=mode)
    assert isinstance(sky_region, SkyRegion)


@pytest.mark.parametrize('region', PIXEL_REGIONS)
def test_pix_to_shapely(region):
    from shapely.geometry import BaseGeometry
    shape = region.to_shapely()
    assert isinstance(shape, BaseGeometry)


@pytest.mark.parametrize(('region', 'mode'), itertools.product(PIXEL_REGIONS, MASK_MODES))
def test_pix_to_mask(region, mode):
    mask = region.to_mask(mode=mode)
    assert isinstance(mask, np.ndarray)
    assert mask.ndim == 2


@pytest.mark.parametrize('region', SKY_REGIONS)
def test_sky_in(region):
    ICRS(1 * u.deg,1 * u.deg) in region


@pytest.mark.parametrize('region', SKY_REGIONS)
def test_sky_area(region):
    area = region.area
    assert isinstance(area, u.Quantity)


@pytest.mark.parametrize(('region', 'mode'), itertools.product(SKY_REGIONS, SKYPIX_MODES))
def test_sky_to_pix(region, mode):
    wcs = WCS(naxis=2)
    pix_region = region.to_pixel(mode=mode, wcs=wcs)
    assert isinstance(pix_region, PixelRegion)
