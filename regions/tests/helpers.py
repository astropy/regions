# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.modeling import models
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions.core import PixCoord, Region

WCS_CENTER = SkyCoord(100 * u.deg, 30 * u.deg)
WCS_CDELT_ARCSEC = 0.1


def make_simple_wcs(skycoord, resolution, size, *, rotation_deg=0):
    """
    Create a simple TAN WCS object.

    Parameters
    ----------
    skycoord : `~astropy.coordinates.SkyCoord`
        The sky coordinate for the center of the WCS.

    resolution : `~astropy.units.Quantity`
        The pixel scale as an angular Quantity.

    size : int
        The size of the WCS in pixels (assumed square).

    rotation_deg : float, optional
        Counter-clockwise rotation angle in degrees. Default is 0
        (axis-aligned).

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        A WCS object with the specified parameters.
    """
    crpix = (size + 1) / 2
    cdelt = resolution.to(u.deg).value
    skycoord_icrs = skycoord.transform_to('icrs')
    ra = skycoord_icrs.ra.degree
    dec = skycoord_icrs.dec.degree

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [crpix, crpix]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    if rotation_deg == 0:
        wcs.wcs.cdelt = np.array([-cdelt, cdelt])
    else:
        theta = np.radians(rotation_deg)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        wcs.wcs.cd = [[-cdelt * cos_t, cdelt * sin_t],
                      [cdelt * sin_t, cdelt * cos_t]]

    return wcs


def make_sip_wcs():
    """
    Create a simple TAN WCS with small SIP distortion terms.

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        A WCS object with SIP distortion.
    """
    header = Header()
    header['NAXIS'] = 2
    header['NAXIS1'] = 20
    header['NAXIS2'] = 20
    header['CRPIX1'] = 10.5
    header['CRPIX2'] = 10.5
    header['CRVAL1'] = 100.0
    header['CRVAL2'] = 30.0
    header['CTYPE1'] = 'RA---TAN-SIP'
    header['CTYPE2'] = 'DEC--TAN-SIP'
    cdelt_deg = WCS_CDELT_ARCSEC / 3600.0
    header['CD1_1'] = -cdelt_deg
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = cdelt_deg
    header['A_ORDER'] = 2
    header['A_2_0'] = 1e-6
    header['B_ORDER'] = 2
    header['B_0_2'] = 1e-6
    return WCS(header)


def make_gwcs():
    """
    Create a simple GWCS object with a TAN projection.

    Returns
    -------
    wcs : `~gwcs.WCS`
        A GWCS object with a TAN projection centred at RA=100, Dec=30.
    """
    from gwcs import WCS as GWCS
    from gwcs import coordinate_frames as cf

    cdelt_deg = WCS_CDELT_ARCSEC / 3600.0
    shift_x = models.Shift(-10.5)
    shift_y = models.Shift(-10.5)
    scale = models.Scale(cdelt_deg) & models.Scale(cdelt_deg)
    rotation = models.AffineTransformation2D(
        matrix=np.array([[-1, 0], [0, 1]]),
        translation=[0, 0],
    )
    tan = models.Pix2Sky_TAN()
    celestial_rotation = models.RotateNative2Celestial(100.0, 30.0, 180.0)

    transform = ((shift_x & shift_y) | scale | rotation
                 | tan | celestial_rotation)

    detector = cf.Frame2D(name='detector', unit=(u.pix, u.pix))
    sky = cf.CelestialFrame(name='sky', axes_order=(0, 1),
                            reference_frame='icrs',
                            axes_names=('lon', 'lat'))

    return GWCS(forward_transform=transform,
                input_frame=detector, output_frame=sky)


def assert_skycoord_allclose(skycoord1, skycoord2, **kwargs):
    """
    Test that two SkyCoord objects are nearly equal.

    Parameters
    ----------
    skycoord1 : `~astropy.coordinates.SkyCoord`
        First SkyCoord object to compare.

    skycoord2 : `~astropy.coordinates.SkyCoord`
        Second SkyCoord object to compare.

    **kwargs : dict
        Additional keyword arguments passed to
        `astropy.tests.assert_quantity_allclose`.
    """
    for attr in (skycoord1._extra_frameattr_names
                 | skycoord2._extra_frameattr_names):
        if not skycoord1.frame._frameattr_equiv(getattr(skycoord1, attr),
                                                getattr(skycoord2, attr)):
            raise ValueError('Extra frame attributes are not equivalent')

    assert skycoord1.is_equivalent_frame(skycoord2)

    for attr in (set(skycoord1._sky_coord_frame._data.components)
                 | set(skycoord2._sky_coord_frame._data.components)):
        sdata1 = skycoord1._sky_coord_frame._data
        sdata2 = skycoord2._sky_coord_frame._data
        assert_quantity_allclose(getattr(sdata1, attr), getattr(sdata2, attr),
                                 **kwargs)


def assert_region_allclose(region1, region2, **kwargs):
    """
    Test that two Region objects have parameters which are nearly equal.

    Meta and Visual properties are matched identically.

    Parameters
    ----------
    region1 : `~regions.core.Region`
        First Region object to compare.

    region2 : `~regions.core.Region`
        Second Region object to compare.

    **kwargs : dict
        Additional keyword arguments passed to
        `astropy.tests.assert_quantity_allclose` or similar functions
        for comparing specific parameter types.
    """
    if not (isinstance(region1, Region) and isinstance(region2, Region)):
        raise TypeError('Both inputs must be Region instances')

    # check that both have identical parameters
    region1_params = list(region1._params)
    region2_params = list(region2._params)
    if region1_params != region2_params:
        raise ValueError('Inputs do not have the same parameters')

    assert region1.meta == region2.meta
    assert region1.visual == region2.visual

    # now check the parameter values
    # Note that Quantity comparisons allow for different units
    # if they directly convertible (e.g., 1. * u.deg == 60. * u.arcmin)
    for param in region1_params:
        value1 = getattr(region1, param)
        value2 = getattr(region2, param)
        if isinstance(value1, SkyCoord):
            assert_skycoord_allclose(value1, value2, **kwargs)
        elif isinstance(value1, PixCoord):
            assert_allclose(value1.xy, value2.xy)
        elif isinstance(value1, u.Quantity):
            assert_quantity_allclose(value1, value2, **kwargs)
        elif isinstance(value1, str):
            assert value1 == value2
        else:
            assert_allclose(value1, value2)
