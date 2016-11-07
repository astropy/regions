# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import pytest, assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.io.fits import getheader
from astropy.wcs import WCS
from ...core import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion

try:
    import matplotlib

    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False

try:
    import wcsaxes

    HAS_WCSAXES = True
except:
    HAS_WCSAXES = False


def test_circle_pixel():
    center = PixCoord(3, 4)
    reg = CirclePixelRegion(center, 2)

    assert str(reg) == 'CirclePixelRegion\ncenter: PixCoord(x=3, y=4)\nradius: 2'


def test_circle_sky():
    center = SkyCoord(3 * u.deg, 4 * u.deg)
    reg = CircleSkyRegion(center, 2 * u.arcsec)

    assert str(reg) == 'CircleSkyRegion\ncenter: <SkyCoord (ICRS): (ra, dec) in deg\n    (3.0, 4.0)>\nradius: 2.0 arcsec'


def test_transformation():
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
    skycircle = CircleSkyRegion(skycoord, 2 * u.arcsec)

    headerfile = get_pkg_data_filename('data/example_header.fits')
    h = getheader(headerfile)
    wcs = WCS(h)

    pixcircle = skycircle.to_pixel(wcs)

    assert_allclose(pixcircle.center.x, -50.5)
    assert_allclose(pixcircle.center.y, 299.5)
    assert_allclose(pixcircle.radius, 0.027777777777828305)

    skycircle2 = pixcircle.to_sky(wcs)

    assert_quantity_allclose(skycircle.center.data.lon, skycircle2.center.data.lon)
    assert_quantity_allclose(skycircle.center.data.lat, skycircle2.center.data.lat)
    assert_quantity_allclose(skycircle2.radius, skycircle.radius)


# Todo : restructure test to use same circle everywhere
# see https://github.com/astropy/regions/pull/20
@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.skipif('not HAS_WCSAXES')
def test_plot():
    import matplotlib.pyplot as plt

    skycoord = SkyCoord(1.983333 * u.deg, -1.99 * u.deg, frame='galactic')
    c = CircleSkyRegion(skycoord, 20 * u.arcsec)
    c.visual.update(color='red')
    headerfile = get_pkg_data_filename('data/example_header.fits')
    h = getheader(headerfile)
    wcs = WCS(h)
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
    p = c.as_patch(ax, alpha=0.6)

    assert_allclose(p.center[0], skycoord.icrs.ra.value)
    assert_allclose(p.center[1], skycoord.icrs.dec.value)
    assert p.get_facecolor() == (1, 0, 0, 0.6)

    c.plot(ax, alpha=0.6)



def test_to_mask():
    center = PixCoord(3, 4)
    reg = CirclePixelRegion(center, 2)
    mask = reg.to_mask(mode='exact')
    assert_allclose(np.sum(mask.data), 4 * np.pi)
