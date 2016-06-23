# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import astropy.units as u
from astropy.coordinates import SkyCoord
from ...shapes import CircleSkyRegion


def test_compound():
    skycoord1 = SkyCoord(0 * u.deg, 0 * u.deg, frame='galactic')
    c1 = CircleSkyRegion(skycoord1, 1 * u.deg)

    skycoord2 = SkyCoord(1 * u.deg, 1 * u.deg, frame='galactic')
    c2 = CircleSkyRegion(skycoord2, 0.5 * u.deg)

    test_coord1 = SkyCoord(1.2 * u.deg, 1.2 * u.deg, frame='galactic')
    test_coord2 = SkyCoord(0.5 * u.deg, 0.5 * u.deg, frame='galactic')
    test_coord3 = SkyCoord(0.7 * u.deg, 0.7 * u.deg, frame='galactic')
    test_coord4 = SkyCoord(2 * u.deg, 5 * u.deg, frame='galactic')

    assert test_coord1 in c2 and test_coord1 not in c1
    assert test_coord2 not in c2 and test_coord2 in c1
    assert test_coord3 in c1 and test_coord3 in c2
    assert test_coord4 not in c2 and test_coord4 not in c1

    coords = SkyCoord([test_coord1, test_coord2, test_coord3, test_coord4], frame='galactic')

    union = c1 | c2
    assert (union.contains(coords) == [True, True, True, False]).all()

    intersection = c1 & c2
    assert (intersection.contains(coords) == [False, False, True, False]).all()

    diff = c1 ^ c2
    assert (diff.contains(coords) == [True, True, False, False]).all()

    c3 = CircleSkyRegion(test_coord4, 0.1 * u.deg)

    union = c1 | c2 | c3
    assert (union.contains(coords) == [True, True, True, True]).all()

    intersection = c1 & c2 & c3
    assert (intersection.contains(coords) == [False, False, False, False]).all()

    diff = c1 ^ c2 ^ c3
    assert (diff.contains(coords) == [True, True, False, True]).all()
