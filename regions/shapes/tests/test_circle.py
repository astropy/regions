from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_filename
from astropy.io.fits import getheader
from astropy.wcs import WCS

from ..circle import CirclePixelRegion, CircleSkyRegion
from ...core import PixCoord


def test_init_pixel():
    pixcoord = PixCoord(3, 4)
    c = CirclePixelRegion(pixcoord, 2)

def test_init_sky():
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    c = CircleSkyRegion(skycoord, 2 * u.arcsec)

# Todo : restructure test to use same circle everywhere
def test_plot():
    import matplotlib.pyplot as plt

    pixcoord = PixCoord(3, 4)
    c = CirclePixelRegion(pixcoord, 2)
    headerfile = get_pkg_data_filename('example_header.fits')
    h = getheader(headerfile)
    wcs = WCS(h)
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs)
    c.to_mpl_patch(ax=ax)