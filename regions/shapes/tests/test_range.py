# Licensed under a 3-clause BSD style license - see LICENSE.rst
from numpy.testing import assert_equal
import astropy.units as u
from astropy.coordinates import SkyCoord
from ..range import RangeSphericalRegion

# TODO: use test cases with finite & infinite values from here:
#  http://www.ivoa.net/documents/SIA/20151223/REC-SIA-2.0-20151223.html#toc12

def test_range_default():
    region = RangeSphericalRegion()
    assert region.contains(SkyCoord("0d 0d"))

def test_range_finite():
    region = RangeSphericalRegion(12 * u.deg, 12.5 * u.deg, 34 * u.deg, 36.0 * u.deg)
    assert region.contains(SkyCoord("12d 34d"))
    assert not region.contains(SkyCoord("11d 34d"))
    assert not region.contains(SkyCoord("13d 34d"))
    assert not region.contains(SkyCoord("12d 33d"))
    assert not region.contains(SkyCoord("12d 37d"))

    mask = region.contains(SkyCoord([12, 13] * u.deg, [34, 34] * u.deg))
    assert_equal(mask, [True, False])
