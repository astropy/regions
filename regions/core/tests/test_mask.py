# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_allclose
import pytest

from ..bounding_box import BoundingBox
from ..mask import RegionMask
from ..pixcoord import PixCoord
from ...shapes import CirclePixelRegion

try:
    import matplotlib    # noqa
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


POSITIONS = [(-20, -20), (-20, 20), (20, -20), (60, 60)]


def test_mask_input_shapes():
    with pytest.raises(ValueError):
        mask_data = np.ones((10, 10))
        bbox = BoundingBox(5, 10, 5, 10)
        RegionMask(mask_data, bbox)


def test_mask_array():
    mask_data = np.ones((10, 10))
    bbox = BoundingBox(5, 15, 5, 15)
    mask = RegionMask(mask_data, bbox)
    data = np.array(mask)
    assert_allclose(data, mask.data)


def test_mask_cutout_shape():
    mask_data = np.ones((10, 10))
    bbox = BoundingBox(5, 15, 5, 15)
    mask = RegionMask(mask_data, bbox)

    with pytest.raises(ValueError):
        mask.cutout(np.arange(10))

    with pytest.raises(ValueError):
        mask._overlap_slices((10,))

    with pytest.raises(ValueError):
        mask.to_image((10,))


def test_mask_cutout_copy():
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(25, 25), radius=10.)
    mask = aper.to_mask()
    cutout = mask.cutout(data, copy=True)
    data[25, 25] = 100.
    assert cutout[10, 10] == 1.


@pytest.mark.parametrize('position', POSITIONS)
def test_mask_cutout_no_overlap(position):
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(position[0], position[1]), radius=10.)
    mask = aper.to_mask()

    cutout = mask.cutout(data)
    assert cutout is None

    weighted_data = mask.multiply(data)
    assert weighted_data is None

    image = mask.to_image(data.shape)
    assert image is None


@pytest.mark.parametrize('position', POSITIONS)
def test_mask_cutout_partial_overlap(position):
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(position[0], position[1]), radius=30.)
    mask = aper.to_mask()

    cutout = mask.cutout(data)
    assert cutout.shape == mask.shape

    weighted_data = mask.multiply(data)
    assert weighted_data.shape == mask.shape

    image = mask.to_image(data.shape)
    assert image.shape == data.shape


def test_mask_nan_in_bbox():
    """
    Regression test that non-finite data values outside of the mask but
    within the bounding box are set to zero.
    """

    data = np.ones((101, 101))
    data[33, 33] = np.nan
    data[67, 67] = np.inf
    data[33, 67] = -np.inf
    data[22, 22] = np.nan
    data[22, 23] = np.inf

    radius = 20.
    reg1 = CirclePixelRegion(PixCoord(50, 50), radius)
    reg2 = CirclePixelRegion(PixCoord(5, 5), radius)

    wdata1 = reg1.to_mask(mode='exact').multiply(data)
    assert_allclose(np.sum(wdata1), np.pi * radius**2)

    wdata2 = reg2.to_mask(mode='exact').multiply(data)
    assert_allclose(np.sum(wdata2), 561.6040111923013)
