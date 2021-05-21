# Licensed under a 3-clause BSD style license - see LICENSE.rst
# (taken from photutils: should probably migrate into astropy.wcs)
from astropy import units as u
import numpy as np

from ..core.pixcoord import PixCoord


def pixel_scale_angle_at_skycoord(skycoord, wcs, offset=1 * u.arcsec):
    """
    Calculate the pixel coordinate and scale and WCS rotation angle at
    the position of a SkyCoord coordinate.

    Parameters
    ----------
    skycoord : `~astropy.coordinates.SkyCoord`
        The SkyCoord coordinate.

    wcs : WCS object
        A world coordinate system (WCS) transformation that
        supports the `astropy shared interface for WCS
        <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ (e.g.,
        `astropy.wcs.WCS`, `gwcs.wcs.WCS`).

    offset : `~astropy.units.Quantity`
        A small angular offset to use to compute the pixel scale and
        position angle.

    Returns
    -------
    pixcoord : `~regions.core.PixCoord`
        The pixel coordinate.

    scale : `~astropy.units.Quantity`
        The pixel scale in arcsec/pixel.

    angle : `~astropy.units.Quantity`
        The angle (in degrees) measured counterclockwise from the
        positive x axis to the "North" axis of the celestial coordinate
        system.

    Notes
    -----
    If distortions are present in the image, the x and y pixel scales
    likely differ.  This function computes a single pixel scale along
    the North/South axis.
    """
    # Convert to pixel coordinates
    x, y = wcs.world_to_pixel(skycoord)
    pixcoord = PixCoord(x=x, y=y)

    # We take a point directly North (i.e., latitude offset) the
    # input sky coordinate and convert it to pixel coordinates,
    # then we use the pixel deltas between the input and offset sky
    # coordinate to calculate the pixel scale and angle.
    skycoord_offset = skycoord.directional_offset_by(0.0, offset)
    x_offset, y_offset = wcs.world_to_pixel(skycoord_offset)

    dx = x_offset - x
    dy = y_offset - y
    scale = offset.to(u.arcsec) / (np.hypot(dx, dy) * u.pixel)
    angle = (np.arctan2(dy, dx) * u.radian).to(u.deg)

    return pixcoord, scale, angle


def assert_angle_or_pixel(name, q):
    """
    Check that ``q`` is either an angular or a pixel
    `~astropy.units.Quantity`.
    """
    if isinstance(q, u.Quantity):
        if q.unit.physical_type == 'angle' or q.unit is u.pixel:
            pass
        else:
            raise ValueError(f"{name} should have angular or pixel units")
    else:
        raise TypeError(f"{name} should be a Quantity instance")


def assert_angle(name, q):
    """
    Check that ``q`` is an angular `~astropy.units.Quantity`.
    """
    if isinstance(q, u.Quantity):
        if q.unit.physical_type == 'angle':
            pass
        else:
            raise ValueError(f"{name} should have angular units")
    else:
        raise TypeError(f"{name} should be a Quantity instance")
