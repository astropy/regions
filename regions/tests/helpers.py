from __future__ import print_function

import numpy as np

from astropy.wcs import WCS
from astropy import units as u


def make_simple_wcs(skycoord, resolution, size):

    crpix = (size + 1) / 2
    cdelt = resolution.to(u.deg).value
    skycoord_icrs = skycoord.transform_to('icrs')
    ra = skycoord_icrs.ra.degree
    dec = skycoord_icrs.dec.degree

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [crpix, crpix]
    wcs.wcs.cdelt = np.array([-cdelt, cdelt])
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    return wcs
