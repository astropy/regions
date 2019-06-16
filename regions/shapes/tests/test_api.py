# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The tests in this file simply check what functionality is currently
implemented and doesn't check anything about correctness.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import itertools
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from ...core.mask import RegionMask
from ...core.core import Region, SkyRegion, PixelRegion
from ...core.pixcoord import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion
from ..ellipse import EllipsePixelRegion, EllipseSkyRegion
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from ..rectangle import RectanglePixelRegion, RectangleSkyRegion
from ..point import PointPixelRegion, PointSkyRegion
from ..annulus import (CircleAnnulusPixelRegion, CircleAnnulusSkyRegion,
                       RectangleAnnulusPixelRegion, RectangleAnnulusSkyRegion,
                       EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion)

PIXEL_REGIONS = [
    CirclePixelRegion(PixCoord(3, 4), radius=5),
    CircleAnnulusPixelRegion(PixCoord(3, 4), 5, 7),
    EllipsePixelRegion(PixCoord(3, 4), width=7, height=5, angle=3 * u.deg),
    EllipseAnnulusPixelRegion(PixCoord(3, 4), 7, 9, 8, 10, angle=3 * u.deg),
    PolygonPixelRegion(PixCoord([1, 4, 3], [2, 4, 4])),
    RectanglePixelRegion(PixCoord(6, 5), width=3, height=5),
    RectangleAnnulusPixelRegion(PixCoord(6, 5), 3, 6, 5, 7),
    PointPixelRegion(PixCoord(1, 2)),
]

SKY_REGIONS = [
    CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), radius=5 * u.deg),
    CircleAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 5 * u.deg, 7 * u.deg),
    EllipseSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), width=7 * u.deg,
                     height=5 * u.deg, angle=3 * u.deg),
    EllipseAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 7 * u.deg, 9 * u.deg,
                            5 * u.deg, 9 * u.deg, angle=3 * u.deg),
    PolygonSkyRegion(SkyCoord([1, 4, 3] * u.deg, [2, 4, 4] * u.deg)),
    RectangleSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg), width=3 * u.deg, height=5 * u.deg),
    RectangleAnnulusSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg), 3 * u.deg, 5 * u.deg,
                              5 * u.deg, 7 * u.deg),
    PointSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg)),
]

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
    # TODO: remove the pytest.skip once polygon area is implemented
    if isinstance(region, PolygonPixelRegion):
        pytest.skip()

    area = region.area
    assert not isinstance(area, u.Quantity)


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_to_sky(region):
    try:
        sky_region = region.to_sky(COMMON_WCS)
        assert isinstance(sky_region, SkyRegion)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize(('region', 'mode'),
                         itertools.product(PIXEL_REGIONS, MASK_MODES),
                         ids=ids_func)
def test_pix_to_mask(region, mode):
    try:
        mask = region.to_mask(mode=mode)
        assert isinstance(mask, RegionMask)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in(region):
    region.contains(SkyCoord(1 * u.deg, 1 * u.deg, frame='icrs'), COMMON_WCS)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in_array(region):
    region.contains(SkyCoord([1, 2, 3] * u.deg, [3, 2, 1] * u.deg, frame='icrs'), COMMON_WCS)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_to_pix(region):
    pix_region = region.to_pixel(wcs=COMMON_WCS)
    assert isinstance(pix_region, PixelRegion)


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_attribute_validation_pixel_regions(region):
    invalid_values = dict(center=[PixCoord([1, 2], [2, 3]), 1,
                                  SkyCoord(1 * u.deg, 2 * u.deg)],
                          radius=[u.Quantity("1deg"), [1], PixCoord(1, 2)],
                          angle=[u.Quantity([1 * u.deg, 2 * u.deg]), 2,
                                 PixCoord(1, 2)],
                          vertices=[u.Quantity("1"), 2, PixCoord(1, 2),
                                    PixCoord([[1, 2]], [[2, 3]])]
                          )
    invalid_values['width'] = invalid_values['radius']
    invalid_values['height'] = invalid_values['radius']
    invalid_values['inner_height'] = invalid_values['radius']
    invalid_values['inner_width'] = invalid_values['radius']
    invalid_values['outer_height'] = invalid_values['radius']
    invalid_values['outer_width'] = invalid_values['radius']
    invalid_values['inner_radius'] = invalid_values['radius']
    invalid_values['outer_radius'] = invalid_values['radius']
    invalid_values['start'] = invalid_values['center']
    invalid_values['end'] = invalid_values['radius']

    for attr in invalid_values:
        if hasattr(region, attr):
            for val in invalid_values.get(attr, None):
                with pytest.raises(ValueError) as err:
                    setattr(region, attr, val)
                assert 'The {} must be'.format(attr) in str(err)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_attribute_validation_sky_regions(region):
    invalid_values = dict(center=[PixCoord([1, 2], [2, 3]), 1,
                                  SkyCoord([1 * u.deg], [2 * u.deg])],
                          radius=[u.Quantity([1 * u.deg, 5 * u.deg]),
                                  [1], SkyCoord(1 * u.deg, 2 * u.deg), 1],
                          angle=[u.Quantity([1 * u.deg, 2 * u.deg]), 2,
                                 SkyCoord(1 * u.deg, 2 * u.deg)],
                          vertices=[u.Quantity("1deg"), 2, SkyCoord(1 * u.deg, 2 * u.deg),
                                    SkyCoord([[1 * u.deg, 2 * u.deg]], [[2 * u.deg, 3 * u.deg]])]
                          )
    invalid_values['width'] = invalid_values['radius']
    invalid_values['height'] = invalid_values['radius']
    invalid_values['inner_height'] = invalid_values['radius']
    invalid_values['inner_width'] = invalid_values['radius']
    invalid_values['outer_height'] = invalid_values['radius']
    invalid_values['outer_width'] = invalid_values['radius']
    invalid_values['inner_radius'] = invalid_values['radius']
    invalid_values['outer_radius'] = invalid_values['radius']
    invalid_values['start'] = invalid_values['center']
    invalid_values['end'] = invalid_values['radius']

    for attr in invalid_values:
        if hasattr(region, attr):
            for val in invalid_values.get(attr, None):
                with pytest.raises(ValueError) as err:
                    setattr(region, attr, val)
                assert 'The {} must be'.format(attr) in str(err)
