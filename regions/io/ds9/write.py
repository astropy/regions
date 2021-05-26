# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..core import to_shape_list

__all__ = ['write_ds9', 'ds9_objects_to_string']


def ds9_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Convert a list of `~regions.Region` objects to a DS9 region string.

    See :ref:`gs-ds9`

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
        A DS9 region string.
    """
    shapelist = to_shape_list(regions, coordsys)
    return shapelist.to_ds9(coordsys, fmt, radunit)


def write_ds9(regions, filename, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Convert a list of `~regions.Region` to a DS9 string and write to a
    file.

    See :ref:`gs-ds9`

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
    output = ds9_objects_to_string(regions, coordsys, fmt, radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
