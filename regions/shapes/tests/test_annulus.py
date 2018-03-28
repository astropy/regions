from __future__ import print_function

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

from ...core import PixCoord
from ...tests.helpers import make_simple_wcs
from ..annulus import CircleAnnulusPixelRegion, CircleAnnulusSkyRegion


def test_init_pixel():

    pixcoord = PixCoord(3, 4)
    annulus = CircleAnnulusPixelRegion(pixcoord, 2, 3)
    assert 'inner' in str(annulus)

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)
    skyannulus = annulus.to_sky(wcs=wcs)
    assert isinstance(skyannulus, CircleAnnulusSkyRegion)

def test_init_sky():

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    inner_radius = 20 * u.arcsec
    outer_radius = 30 * u.arcsec

    annulus = CircleAnnulusSkyRegion(skycoord, inner_radius, outer_radius)

    assert_quantity_allclose(annulus.center.ra, skycoord.ra)
    assert_quantity_allclose(annulus.inner_radius, inner_radius)
    assert_quantity_allclose(annulus.outer_radius, outer_radius)

    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)
    assert not annulus.contains(skycoord, wcs)

    test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
    assert not annulus.contains(test_coord, wcs)

    test_coord = SkyCoord(3 * u.deg, 4.007 * u.deg, frame='icrs')
    assert annulus.contains(test_coord, wcs)

    assert 'Annulus' in str(annulus)
    assert 'inner' in str(annulus)

    pixannulus = annulus.to_pixel(wcs=wcs)
    assert isinstance(pixannulus, CircleAnnulusPixelRegion)
