# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the compound module.
"""

import operator

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Latitude, Longitude, SkyCoord
from numpy.testing import assert_allclose

from regions._utils.examples import _make_example_dataset
from regions.core import (CompoundPixelRegion, CompoundSkyRegion,
                          CompoundSphericalSkyRegion, PixCoord,
                          RegionBoundingBox)
from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            CircleSphericalSkyRegion, PolygonPixelRegion,
                            RangeSphericalSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def fixture_wcs():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = _make_example_dataset(config=config)
    return dataset.wcs


class TestCompoundPixel:
    """
    Test compound pixel regions.
    """

    # Two circles that overlap in one column
    c1 = CirclePixelRegion(PixCoord(5, 5), 4)
    c2 = CirclePixelRegion(PixCoord(11, 5), 4)

    def test_copy(self):
        reg = (self.c1 | self.c2).copy()
        assert isinstance(reg, CompoundPixelRegion)
        assert reg.operator is operator.or_

        assert isinstance(reg.region1, CirclePixelRegion)
        assert reg.region1.center.xy == (5, 5)
        assert reg.region1.radius == 4
        assert reg.region1.meta == {}
        assert reg.region1.visual == {}

        assert isinstance(reg.region2, CirclePixelRegion)
        assert reg.region2.center.xy == (11, 5)
        assert reg.region2.radius == 4
        assert reg.region2.meta == {}
        assert reg.region2.visual == {}

    def test_rotate(self):
        reg = (self.c1 | self.c2).rotate(PixCoord(0, 0), 360 * u.deg)
        assert isinstance(reg, CompoundPixelRegion)
        assert reg.operator is operator.or_

        assert isinstance(reg.region1, CirclePixelRegion)
        assert_allclose(reg.region1.center.xy, (5, 5))
        assert reg.region1.radius == 4
        assert reg.region1.meta == {}
        assert reg.region1.visual == {}

        assert isinstance(reg.region2, CirclePixelRegion)
        assert_allclose(reg.region2.center.xy, (11, 5))
        assert reg.region2.radius == 4
        assert reg.region2.meta == {}
        assert reg.region2.visual == {}

    def test_union(self, wcs):
        union = self.c1 | self.c2
        assert isinstance(union, CompoundPixelRegion)
        mask = union.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

        union_to_sphsky = union.to_spherical_sky(wcs)
        assert isinstance(union_to_sphsky, CompoundSphericalSkyRegion)

        union_to_sky = union.to_sky(wcs)
        assert isinstance(union_to_sky, CompoundSkyRegion)

    def test_intersection(self):
        intersection = self.c1 & self.c2
        mask = intersection.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_symdiff(self):
        symdiff = self.c1 ^ self.c2
        mask = symdiff.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 0, 0, 0, 0, 0, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

    def test_to_mask(self):
        # fixes #168
        pixcoord3 = PixCoord(1, 1)
        c3 = CirclePixelRegion(pixcoord3, 4)
        union2 = self.c1 | c3
        mask1 = union2.to_mask()

        pixcoord4 = PixCoord(9, 9)
        c4 = CirclePixelRegion(pixcoord4, 4)
        union3 = self.c1 | c4
        mask2 = union3.to_mask()

        # mask1 and mask2 should be equal
        assert_allclose(mask1.data, mask2.data)

        ref_data = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        assert_allclose(mask1.data, ref_data)

    def test_bounding_box(self):
        bbox = (self.c1 | self.c2).bounding_box
        assert bbox == RegionBoundingBox(1, 16, 1, 10)

    @pytest.mark.parametrize('method', ['contains', 'covers'])
    def test_contains_scalar_polygon(self, method):
        """
        Test that a compound region with a polygon component returns a
        scalar for a scalar position.
        """
        poly = PolygonPixelRegion(PixCoord([8, 14, 14, 8], [2, 2, 8, 8]))
        reg = self.c1 | poly
        result = getattr(reg, method)(PixCoord(12, 5))
        assert np.ndim(result) == 0
        assert result

    def test_invalid(self):
        with pytest.raises(ValueError):
            self.c1 | 1.0
        with pytest.raises(ValueError):
            self.c1 | PixCoord(1.0, 2.0)


def test_compound_sky():
    skycoord1 = SkyCoord(0 * u.deg, 0 * u.deg, frame='galactic')
    c1 = CircleSkyRegion(skycoord1, 1 * u.deg)

    skycoord2 = SkyCoord(1 * u.deg, 1 * u.deg, frame='galactic')
    c2 = CircleSkyRegion(skycoord2, 0.5 * u.deg)

    test_coord1 = SkyCoord(1.2 * u.deg, 1.2 * u.deg, frame='galactic')
    test_coord2 = SkyCoord(0.5 * u.deg, 0.5 * u.deg, frame='galactic')
    test_coord3 = SkyCoord(0.7 * u.deg, 0.7 * u.deg, frame='galactic')
    test_coord4 = SkyCoord(2 * u.deg, 5 * u.deg, frame='galactic')

    wcs = make_simple_wcs(skycoord1, 0.1 * u.deg, 20)

    assert c2.contains(test_coord1, wcs) and not c1.contains(test_coord1, wcs)
    assert not c2.contains(test_coord2, wcs) and c1.contains(test_coord2, wcs)
    assert c1.contains(test_coord3, wcs) and c2.contains(test_coord3, wcs)
    assert (not c2.contains(test_coord4, wcs)
            and not c1.contains(test_coord4, wcs))

    coords = SkyCoord([test_coord1, test_coord2, test_coord3, test_coord4],
                      frame='galactic')

    union = c1 | c2
    assert (union.contains(coords, wcs) == [True, True, True, False]).all()

    intersection = c1 & c2
    assert ((intersection.contains(coords, wcs)
             == [False, False, True, False]).all())

    diff = c1 ^ c2
    assert (diff.contains(coords, wcs) == [True, True, False, False]).all()

    c3 = CircleSkyRegion(test_coord4, 0.1 * u.deg)
    union = c1 | c2 | c3
    assert (union.contains(coords, wcs) == [True, True, True, True]).all()

    intersection = c1 & c2 & c3
    assert ((intersection.contains(coords, wcs)
             == [False, False, False, False]).all())

    diff = c1 ^ c2 ^ c3
    assert (diff.contains(coords, wcs) == [True, True, False, True]).all()
    assert 'Compound' in str(union)

    # Test keeping the input meta
    meta = {'test_meta': True}
    compound = CompoundSkyRegion(c1, c2, operator.and_, meta=meta)
    assert compound.meta['test_meta']

    union_to_sphsky = union.to_spherical_sky(wcs=wcs)
    assert isinstance(union_to_sphsky, CompoundSphericalSkyRegion)

    union_to_pixel = union.to_pixel(wcs)
    assert isinstance(union_to_pixel, CompoundPixelRegion)


class TestCompoundSphericalSky:
    """
    Test compound spherical sky regions.
    """

    # Test regions and coordinates
    skycoord1 = SkyCoord(0 * u.deg, 0 * u.deg, frame='galactic')
    c1 = CircleSphericalSkyRegion(skycoord1, 1 * u.deg)

    skycoord2 = SkyCoord(1 * u.deg, 1 * u.deg, frame='galactic')
    c2 = CircleSphericalSkyRegion(skycoord2, 0.5 * u.deg)

    test_coord1 = SkyCoord(1.2 * u.deg, 1.2 * u.deg, frame='galactic')
    test_coord2 = SkyCoord(0.5 * u.deg, 0.5 * u.deg, frame='galactic')
    test_coord3 = SkyCoord(0.7 * u.deg, 0.7 * u.deg, frame='galactic')
    test_coord4 = SkyCoord(2 * u.deg, 5 * u.deg, frame='galactic')

    coords = SkyCoord([test_coord1, test_coord2, test_coord3, test_coord4],
                      frame='galactic')

    c3 = CircleSphericalSkyRegion(test_coord4, 0.1 * u.deg)

    skycoord4 = SkyCoord(1 * u.deg, 89 * u.deg, frame='galactic')
    c4 = CircleSphericalSkyRegion(skycoord4, 2.5 * u.deg)

    def test_contains_components(self):
        # Check test point positions relative to component regions
        assert (self.c2.contains(self.test_coord1)
                and not self.c1.contains(self.test_coord1))
        assert (not self.c2.contains(self.test_coord2)
                and self.c1.contains(self.test_coord2))
        assert (self.c1.contains(self.test_coord3)
                and self.c2.contains(self.test_coord3))
        assert (not self.c2.contains(self.test_coord4)
                and not self.c1.contains(self.test_coord4))

    def test_union(self):
        union = self.c1 | self.c2
        assert (union.contains(self.coords) == [True, True, True, False]).all()

        union2 = self.c1 | self.c2 | self.c3
        assert (union2.contains(self.coords) == [True, True, True, True]).all()

        assert 'Compound' in str(union)

    def test_intersection(self):
        intersection = self.c1 & self.c2
        assert ((intersection.contains(self.coords)
                == [False, False, True, False]).all())

        intersection2 = self.c1 & self.c2 & self.c3
        assert ((intersection2.contains(self.coords)
                == [False, False, False, False]).all())

    def test_symdiff(self):
        diff = self.c1 ^ self.c2
        assert (diff.contains(self.coords) == [True, True, False, False]).all()

        diff2 = self.c1 ^ self.c2 ^ self.c3
        assert (diff2.contains(self.coords) == [True, True, False, True]).all()

    def test_transformation(self, wcs):
        union = self.c1 | self.c2

        union_to_sky = union.to_sky(wcs=wcs)
        assert isinstance(union_to_sky, CompoundSkyRegion)

        union_to_pixel = union.to_pixel(wcs)
        assert isinstance(union_to_pixel, CompoundPixelRegion)

    def test_frame_transformation(self):
        union = self.c1 | self.c2
        assert isinstance(union.transform_to('icrs'),
                          CompoundSphericalSkyRegion)

    def test_discretize_boundary(self):
        union = self.c1 | self.c2

        assert isinstance(union.discretize_boundary(n_vertices=2),
                          CompoundSphericalSkyRegion)

    def test_bounding_circle(self):
        union = self.c1 | self.c2
        assert isinstance(union.bounding_circle, CircleSphericalSkyRegion)

    def test_bounding_lonlat(self):
        union = self.c1 | self.c2

        bound_lonlat = union.bounding_lonlat
        assert isinstance(bound_lonlat[0], Longitude)
        assert isinstance(bound_lonlat[1], Latitude)

        # Test no longitude bounds on union
        union3 = self.c1 | self.c4
        bound_lonlat3 = union3.bounding_lonlat
        assert bound_lonlat3[0] is None
        assert isinstance(bound_lonlat3[1], Latitude)

    def test_compound_range(self):
        """
        Regression test: confirm a compound region can be created
        using RangeSphericalSkyRegion instances.
        """
        reg1 = RangeSphericalSkyRegion(
            frame='icrs',
            latitude_range=[0 * u.deg, 1 * u.deg],
        )
        reg2 = RangeSphericalSkyRegion(
            frame='icrs',
            latitude_range=[0.5 * u.deg, 1.5 * u.deg],
        )
        union = reg1 & reg2

        assert isinstance(union, CompoundSphericalSkyRegion)
        assert union.frame.name == 'icrs'
