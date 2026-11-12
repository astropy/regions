# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file sets up detailed tests for computing masks.

The masks are checked against independent calculations. An exact
mask must sum to the analytic area of the region, and each of its
pixels must match the fraction of a fine grid of points within that
pixel that the region contains. A center or subpixel mask must equal
the fraction of the pixel or subpixel centers that the region
contains.
"""
import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from regions.core import PixCoord, RegionMask
from regions.shapes.circle import CirclePixelRegion
from regions.shapes.ellipse import EllipsePixelRegion
from regions.shapes.polygon import PolygonPixelRegion
from regions.shapes.rectangle import RectanglePixelRegion

POLY_X = np.array([-2.334, 3.631, 1.122, -1.341])
POLY_Y = np.array([-3.121, -2.118, 2.987, 1.759])

REGIONS = {
    'circ': CirclePixelRegion(PixCoord(3.981987, 4.131378), radius=3.3411),
    'elli': EllipsePixelRegion(PixCoord(5.981987, 4.131378), width=10.4466,
                               height=6.6822, angle=12 * u.deg),
    'rect': RectanglePixelRegion(PixCoord(5.981987, 4.131378), width=7.2233,
                                 height=4.3411, angle=12 * u.deg),
    'poly': PolygonPixelRegion(PixCoord(POLY_X, POLY_Y)),
}

# The analytic area of each region. The polygon area is computed with
# the shoelace formula.
AREAS = {
    'circ': np.pi * 3.3411**2,
    'elli': np.pi * 10.4466 * 6.6822 / 4,
    'rect': 7.2233 * 4.3411,
    'poly': 0.5 * abs(np.dot(POLY_X, np.roll(POLY_Y, 1))
                      - np.dot(POLY_Y, np.roll(POLY_X, 1))),
}


def sampled_fraction(region, bbox, nsamples):
    """
    Return the fraction of an ``nsamples`` x ``nsamples`` grid of
    points within each pixel of ``bbox`` that the region contains.

    The points are the centers of the subpixels, so ``nsamples=1``
    samples only the pixel centers.
    """
    offsets = (np.arange(nsamples) + 0.5) / nsamples - 0.5
    xpix = bbox.ixmin + np.arange(bbox.shape[1])
    ypix = bbox.iymin + np.arange(bbox.shape[0])
    xsamples = (xpix[:, np.newaxis] + offsets).ravel()
    ysamples = (ypix[:, np.newaxis] + offsets).ravel()
    xgrid, ygrid = np.meshgrid(xsamples, ysamples)
    inside = region.contains(PixCoord(xgrid, ygrid)).astype(float)
    inside = inside.reshape(bbox.shape[0], nsamples,
                            bbox.shape[1], nsamples)
    return inside.mean(axis=(1, 3))


@pytest.mark.parametrize('name', REGIONS)
def test_to_mask_exact(name):
    region = REGIONS[name]
    mask = region.to_mask(mode='exact')
    assert isinstance(mask, RegionMask)
    assert_allclose(mask.data.sum(), AREAS[name], rtol=1e-12)
    assert np.all(mask.data >= 0)
    assert np.all(mask.data <= 1 + 1e-12)

    # With 200 x 200 points per pixel, the sampled fraction of a
    # partially covered pixel differs from its exact overlap by less
    # than about 6e-4 for these regions.
    assert_allclose(mask.data, sampled_fraction(region, mask.bbox, 200),
                    rtol=0, atol=1e-3)


@pytest.mark.parametrize('name', REGIONS)
def test_to_mask_center(name):
    region = REGIONS[name]
    mask = region.to_mask(mode='center')
    assert isinstance(mask, RegionMask)
    assert_equal(mask.data, sampled_fraction(region, mask.bbox, 1))


@pytest.mark.parametrize('subpixels', [1, 5, 18])
@pytest.mark.parametrize('name', REGIONS)
def test_to_mask_subpixels(name, subpixels):
    region = REGIONS[name]
    mask = region.to_mask(mode='subpixels', subpixels=subpixels)
    assert isinstance(mask, RegionMask)
    assert_equal(mask.data,
                 sampled_fraction(region, mask.bbox, subpixels))
