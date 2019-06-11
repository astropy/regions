# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from ..core import to_shape_list

__all__ = [
    'write_crtf',
    'crtf_objects_to_string',
]


def crtf_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Converts a `list` of `~regions.Region` to CRTF region string.

    Parameters
    ----------
    regions : `list`
        List of `~regions.Region` objects
    coordsys : `str`, optional
        Astropy Coordinate system that overrides the coordinate system frame for
        all regions. Default is 'fk5'.
    fmt : `str`, optional
        A python string format defining the output precision. Default is .6f,
        which is accurate to 0.0036 arcseconds.
    radunit : `str`, optional
        This denotes the unit of the radius. Default is deg (degrees)

    Returns
    -------
    region_string : `str`
        CRTF region string

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from regions import CircleSkyRegion, crtf_objects_to_string
    >>> reg_sky = CircleSkyRegion(SkyCoord(1 * u.deg, 2 * u.deg), 5 * u.deg)
    >>> print(crtf_objects_to_string([reg_sky]))
    #CRTF
    circle[[1.000007deg, 2.000002deg], 5.000000deg], coord=FK5,

    """
    shapelist = to_shape_list(regions, coordsys)
    return shapelist.to_crtf(coordsys, fmt, radunit)


def write_crtf(regions, filename, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Converts a `list` of `~regions.Region` to CRTF string and write to file.

    See :ref:`gs-crtf`

    Parameters
    ----------
    regions : `list`
        List of `~regions.Region` objects
    filename : `str`
        Filename in which the string is to be written. Default is 'new.crtf'
    coordsys : `str`, optional
        Astropy Coordinate system that overrides the coordinate frames of all
        regions. Default is 'fk5'.
    fmt : `str`, optional
        A python string format defining the output precision. Default is .6f,
        which is accurate to 0.0036 arcseconds.
    radunit : `str`, optional
        This denotes the unit of the radius. Default is deg (degrees)
    """
    output = crtf_objects_to_string(regions, coordsys, fmt, radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
