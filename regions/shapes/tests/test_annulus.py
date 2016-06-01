from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import pytest
from astropy.utils.data import get_pkg_data_filename
from astropy.io.fits import getheader
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from ..annulus import AnnulusPixelRegion, AnnulusSkyRegion
from ...core import PixCoord

def test_init_pixel():
    pixcoord = PixCoord(3, 4)
    c = AnnulusPixelRegion(pixcoord, 2, 3)

def test_init_sky():
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    c = AnnulusSkyRegion(skycoord, 2 * u.arcsec, 3 * u.arcsec)
