# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..core import to_shape_list

__all__ = ['write_crtf', 'crtf_objects_to_string']


def crtf_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Convert a list of `~regions.Region` objects to a CRTF region string.

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    coordsys : str, optional
        An Astropy coordinate system that overrides the coordinate
        system frame for all regions.

    fmt : str, optional
        A python string format defining the output precision. Default is
        '.6f', which is accurate to 0.0036 arcseconds.

    radunit : str, optional
        The unit of the radius.

    Returns
    -------
    region_string : str
        A CRTF region string.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from regions import CircleSkyRegion, crtf_objects_to_string
    >>> reg_sky = CircleSkyRegion(SkyCoord(1 * u.deg, 2 * u.deg), 5 * u.deg)
    >>> print(crtf_objects_to_string([reg_sky]))
    #CRTFv0
    global coord=J2000
    circle[[1.000007deg, 2.000002deg], 5.000000deg]
    """
    shapelist = to_shape_list(regions, coordsys)
    return shapelist.to_crtf(coordsys, fmt, radunit)


def write_crtf(regions, filename, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Convert a list of `~regions.Region` to a CRTF string and write to a
    file.

    See :ref:`gs-crtf`

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    filename : str
        The filename in which the string is to be written.

    coordsys : str, optional
        An Astropy coordinate system that overrides the coordinate
        frames of all regions.

    fmt : str, optional
        A python string format defining the output precision. Default is
        '.6f', which is accurate to 0.0036 arcseconds.

    radunit : str, optional
        The unit of the radius.
    """
    output = crtf_objects_to_string(regions, coordsys, fmt, radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
