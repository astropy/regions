import math
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from ...core import PixCoord
from ..circle import CirclePixelRegion, CircleSkyRegion

def test_init_pixel():
    pixcoord = PixCoord(3, 4)
    c = CirclePixelRegion(pixcoord, 2)

def test_init_sky():
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg)
    c = CircleSkyRegion(skycoord, 2 * u.arcsec)
