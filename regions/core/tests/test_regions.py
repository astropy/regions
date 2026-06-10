# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the regions module.
"""

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from regions import CirclePixelRegion, CircleSkyRegion, RectanglePixelRegion
from regions.core import PixCoord, Region, Regions


def test_regions_inputs():
    regs = []
    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regs.append(CirclePixelRegion(center, radius=radius))

    regions = Regions(regs)
    assert len(regions) == 4


def test_regions_no_input():
    regions = Regions()
    assert len(regions) == 0

    regions = Regions([])
    assert len(regions) == 0


def test_regions_invalid_input():
    match = "'int' object is not iterable"
    with pytest.raises(TypeError, match=match):
        Regions(1)

    match = 'Input regions must be a list of Region objects'
    with pytest.raises(TypeError, match=match):
        Regions([1])
    with pytest.raises(TypeError, match=match):
        Regions([1, 2, 3])


def test_regions_append():
    regions = Regions()

    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regions.append(CirclePixelRegion(center, radius=radius))
    assert len(regions) == 4

    match = 'Input region must be a Region object'
    with pytest.raises(TypeError, match=match):
        regions.append(1)


def test_regions_extend():
    regions = Regions()

    regs = []
    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regs.append(CirclePixelRegion(center, radius=radius))

    regions.extend(regs)
    assert len(regions) == 4

    regions.extend(regions)
    assert len(regions) == 8

    match = 'Input regions must be a list of Region objects'
    with pytest.raises(TypeError, match=match):
        regions.extend([1])


def test_regions_methods():
    regions = Regions()

    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regions.append(CirclePixelRegion(center, radius=radius))

    reg_slc = regions[0:2]
    assert len(reg_slc) == 2

    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions.insert(0, reg)
    assert len(regions) == 5
    assert regions[0] == reg

    regions2 = regions.copy()
    regions.reverse()
    assert regions.regions == regions2.regions[::-1]

    outreg = regions.pop(-1)
    assert outreg == reg


def test_regions_get_formats():
    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions = Regions([reg])
    tbl = regions.get_formats()
    assert isinstance(tbl, Table)
    assert len(tbl) == 3


def test_regions_repr():
    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions = Regions([reg])
    reg_repr = '<CirclePixelRegion(center=PixCoord(x=0, y=0), radius=1)>'
    regions_str = reg_repr
    regions_repr = f'<Regions([\n  {reg_repr}\n])>'
    assert str(regions) == regions_str
    assert repr(regions) == regions_repr


def test_region_equivalency():
    """
    Test the __eq__ method of Region and its subclasses to ensure that
    it correctly identifies when two regions are equal, and that it uses
    np.allclose for comparing parameters that are numerical in nature.
    """
    reg1 = CirclePixelRegion(PixCoord(0, 0), radius=1)
    reg2 = CirclePixelRegion(PixCoord(0, 0), radius=1)
    assert reg1 == reg2

    # __eq__ uses np.allclose for comparing regions, so this should be True
    reg3 = CirclePixelRegion(PixCoord(0, 0), radius=1 + 1e-10)
    assert reg1 == reg3

    # Different radius should not be equal
    reg_diff = CirclePixelRegion(PixCoord(0, 0), radius=2)
    assert reg1 != reg_diff

    # Test SkyRegion equality
    reg4 = CircleSkyRegion(SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
                           radius=1 * u.arcsec)
    reg5 = CircleSkyRegion(SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
                           radius=1 * u.arcsec)
    assert reg4 == reg5

    reg6 = CircleSkyRegion(SkyCoord(ra=0 * u.deg, dec=0 * u.deg),
                           radius=1 * u.arcsec + 1e-10 * u.arcsec)
    assert reg4 == reg6

    # Different classes should not be equal
    assert reg1 != reg4

    # Test SkyCoord center allclose comparison
    # Small coordinate differences should be considered equal
    reg7 = CircleSkyRegion(SkyCoord(ra=0 * u.deg + 1e-10 * u.deg,
                                    dec=0 * u.deg),
                           radius=1 * u.arcsec)
    assert reg4 == reg7

    # Declination offset
    reg8 = CircleSkyRegion(SkyCoord(ra=0 * u.deg,
                                    dec=0 * u.deg + 1e-10 * u.deg),
                           radius=1 * u.arcsec)
    assert reg4 == reg8

    # Both RA and Dec with small offsets
    reg9 = CircleSkyRegion(SkyCoord(ra=0 * u.deg + 1e-10 * u.deg,
                                    dec=0 * u.deg + 1e-10 * u.deg),
                           radius=1 * u.arcsec)
    assert reg4 == reg9

    # Significantly different SkyCoord should not be equal
    reg10 = CircleSkyRegion(SkyCoord(ra=0.01 * u.deg, dec=0 * u.deg),
                            radius=1 * u.arcsec)
    assert reg4 != reg10

    # SkyCoord with different frames should not be equal
    reg_gal = CircleSkyRegion(SkyCoord(l=0 * u.deg, b=0 * u.deg,
                                       frame='galactic'),
                              radius=1 * u.arcsec)
    assert reg4 != reg_gal


def test_is_close_skycoord_frame_error(monkeypatch):
    """
    Test that if the is_equivalent_frame method of the SkyCoord frame
    raises a ValueError, then _is_close should catch that and return
    False rather than propagating the exception.

    This is to ensure that if the frames of the SkyCoords cannot be
    compared for equivalency, we should not consider the coordinates to
    be close, but we also should not raise an exception.
    """
    def raises_value_error(self, other):
        raise ValueError

    c1 = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
    c2 = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
    monkeypatch.setattr(type(c1.frame), 'is_equivalent_frame',
                        raises_value_error)
    assert not Region._is_close(c1, c2)


def test_is_close_non_comparable():
    """
    Test that if the __eq__ method of the objects being compared
    raises a TypeError, then _is_close should return False rather than
    propagating the exception.
    """
    class BadCompare:
        def __eq__(self, other):
            raise TypeError

    bad = BadCompare()
    assert not Region._is_close(bad, bad)


def test_region_eq_param_mismatch():
    """
    Test that if two regions are of the same class but have mismatched
    _params, then they should not be considered equal.

    This is a sanity check to ensure that the __eq__ method is correctly
    comparing the parameters of the regions.
    """
    reg1 = CirclePixelRegion(PixCoord(0, 0), radius=1)
    reg2 = CirclePixelRegion(PixCoord(0, 0), radius=1)
    reg2._params = ()  # simulate mismatched parameter set
    assert reg1 != reg2


def test_compound_covers():
    """
    Test covers method on compound regions.
    """
    rect = RectanglePixelRegion(center=PixCoord(1, 1), width=2, height=2)
    comp_union = rect.union(rect)
    assert comp_union.covers(PixCoord(1, 1))

    comp_inter = rect.intersection(rect)
    assert comp_inter.covers(PixCoord(2, 2))  # on the boundary
