# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from regions.io.asdf.converters.tests import examples


@pytest.fixture
def test_pixcoord():
    """
    Test that PixCoord can be serialized and deserialized correctly.
    """
    return examples.pixcoord()


@pytest.fixture
def test_pixcoord_array():
    """
    Test that PixCoord can be serialized and deserialized correctly.
    """
    return examples.pixcoord_array()


@pytest.fixture
def test_circle_pixel_region():
    """
    Test that CirclePixelRegion can be serialized and deserialized correctly.
    """
    return examples.circle_pixel_region()


@pytest.fixture
def test_circle_sky_region():
    """
    Test that CircleSkyRegion can be serialized and deserialized correctly.
    """
    return examples.circle_sky_region()


@pytest.fixture
def test_circle_spherical_sky_region():
    """
    Test 'CircleSphericalSkyRegion' can be serialized and deserialized.
    """
    return examples.circle_spherical_sky_region()


@pytest.fixture
def test_meta():
    """
    Test that RegionMeta can be serialized and deserialized correctly.
    """
    return examples.meta()


@pytest.fixture
def test_visual_meta():
    """
    Test that RegionVisual can be serialized and deserialized correctly.
    """
    return examples.visual_meta()


@pytest.fixture
def test_ellipse_pixel_region():
    """
    Test that EllipsePixelRegion can be serialized and deserialized correctly.
    """
    return examples.ellipse_pixel_region()


@pytest.fixture
def test_ellipse_sky_region():
    """
    Test that EllipseSkyRegion can be serialized and deserialized correctly.
    """
    return examples.ellipse_sky_region()
