# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...shapes import CircleSkyRegion, CirclePixelRegion
from ...core import PixCoord, CompoundPixelRegion
from ...tests.helpers import make_simple_wcs
import pytest
from numpy.testing import assert_allclose

def test_compound_pixel():
    # Two circles that overlap in one column
    pixcoord1 = PixCoord(5, 5)
    c1 = CirclePixelRegion(pixcoord1, 4)
    pixcoord2 = PixCoord(11, 5)
    c2 = CirclePixelRegion(pixcoord2, 4)

    union = c1 | c2
    assert isinstance(union, CompoundPixelRegion)
    mask = union.to_mask()
    assert_allclose(mask.data[:,7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
    assert_allclose(mask.data[:,6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

    intersection = c1 & c2
    mask = intersection.to_mask()
    assert_allclose(mask.data[:,7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
    assert_allclose(mask.data[:,6], [0, 0, 0, 0, 0, 0, 0, 0, 0])

    symdiff = c1 ^ c2
    mask = symdiff.to_mask()
    assert_allclose(mask.data[:,7], [0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert_allclose(mask.data[:,6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

    # Circle bigger than map, see #168
    pixcoord3 = PixCoord(1, 1)
    c3 = CirclePixelRegion(pixcoord3, 4)
    union2 = c1 | c3
    with pytest.raises(NotImplementedError):
        mask = union2.to_mask()

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
    assert not c2.contains(test_coord4, wcs) and not c1.contains(test_coord4, wcs)

    coords = SkyCoord([test_coord1, test_coord2, test_coord3, test_coord4], frame='galactic')

    union = c1 | c2
    assert (union.contains(coords, wcs) == [True, True, True, False]).all()

    intersection = c1 & c2
    assert (intersection.contains(coords, wcs) == [False, False, True, False]).all()

    diff = c1 ^ c2
    assert (diff.contains(coords, wcs) == [True, True, False, False]).all()

    c3 = CircleSkyRegion(test_coord4, 0.1 * u.deg)

    union = c1 | c2 | c3
    assert (union.contains(coords, wcs) == [True, True, True, True]).all()

    intersection = c1 & c2 & c3
    assert (intersection.contains(coords, wcs) == [False, False, False, False]).all()

    diff = c1 ^ c2 ^ c3
    assert (diff.contains(coords, wcs) == [True, True, False, True]).all()

    assert 'Compound' in str(union)
