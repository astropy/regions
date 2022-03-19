# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the fits subpackage.
"""

import warnings

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy.utils.data import get_pkg_data_filenames
from astropy.utils.exceptions import AstropyUserWarning
from astropy.table import QTable
from astropy.wcs import WCS
from numpy.testing import assert_allclose
import pytest

from ....core import Regions, PixCoord
from ....shapes import CirclePixelRegion, CircleSkyRegion
from ....tests.helpers import assert_region_allclose
from ..core import FITSParserError


def test_roundtrip(tmpdir):
    filenames = get_pkg_data_filenames('data', pattern='*.fits')

    # AstropyUserWarning will be emitted only for some of the files
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        # Check that all test files are readable
        for filename in filenames:
            regions = Regions.read(filename, format='fits')

    tempfile = tmpdir.join('tmp.fits').strpath
    regions.write(tempfile, format='fits', overwrite=True)
    regions2 = Regions.read(tempfile, format='fits')
    for reg1, reg2 in zip(regions, regions2):
        assert_region_allclose(reg1, reg2)


def test_only_pixel_regions():
    reg_sky = CircleSkyRegion(SkyCoord(1, 2, unit='deg'), 5 * u.deg)
    reg = Regions([reg_sky])

    match = 'Sky regions cannot be serialized'
    with pytest.warns(AstropyUserWarning, match=match):
        result = reg.serialize(format='fits')
        assert len(result) == 0

    reg_pix = CirclePixelRegion(PixCoord(10, 10), 5)
    reg = Regions([reg_sky, reg_pix])
    match = 'Sky regions cannot be serialized'
    with pytest.warns(AstropyUserWarning, match=match):
        result = reg.serialize(format='fits')
        assert len(result) == 1


def test_valid_columns():
    t = QTable([[1, 2, 3]], names=('a'))
    with pytest.raises(FITSParserError) as excinfo:
        Regions.parse(t, format='fits')
        estr = "'a' is not a valid column name"
        assert estr in str(excinfo.value)


def test_valid_row():
    x = [1] * u.pix
    y = [2] * u.pix
    shapes = ['CIRCLE']
    tbl = QTable([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    match = 'Table columns are missing'
    with pytest.warns(AstropyUserWarning, match=match):
        Regions.parse(tbl, format='fits')

    # test invald shape
    shapes = ['INVALID']
    tbl2 = QTable([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    with pytest.raises(FITSParserError) as excinfo:
        Regions.parse(tbl2, format='fits')
        estr = "'invalid' is not a valid FITS region shape"
        assert estr in str(excinfo.value)

    shapes = ['PIE']
    rotang = [[20, 30]] * u.deg
    tbl3 = QTable([x, y, shapes, rotang], names=('X', 'Y', 'SHAPE', 'ROTANG'))
    match = "'pie' is not supported"
    with pytest.warns(AstropyUserWarning, match=match):
        Regions.parse(tbl3, format='fits')
