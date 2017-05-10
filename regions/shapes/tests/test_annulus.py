from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import pytest, assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io.fits import getheader
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from ..annulus import CircleAnnulusPixelRegion, CircleAnnulusSkyRegion
from ...core import PixCoord

def test_init_pixel():
    pixcoord = PixCoord(3, 4)
    annulus = CircleAnnulusPixelRegion(pixcoord, 2, 3)
    assert 'inner' in str(annulus)

def test_init_sky():
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    inner_radius = 20 * u.arcsec
    outer_radius = 30 * u.arcsec
    annulus= CircleAnnulusSkyRegion(skycoord, inner_radius, outer_radius)

    assert_quantity_allclose(annulus.center.ra, skycoord.ra)
    assert_quantity_allclose(annulus.inner_radius, inner_radius)
    assert_quantity_allclose(annulus.outer_radius, outer_radius)

    assert skycoord not in annulus

    test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
    assert test_coord not in annulus

    test_coord = SkyCoord(3 * u.deg, 4.007 * u.deg, frame='icrs')
    assert test_coord in annulus

    assert 'Annulus' in str(annulus)
    assert 'inner' in str(annulus)
