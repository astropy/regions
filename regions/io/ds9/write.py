# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from ..core import to_shape_list

__all__ = [
    'write_ds9',
    'ds9_objects_to_string',
]


def ds9_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Converts list of regions to ds9 region string.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    coordsys : str
        This overrides the coordinate system frame for all regions.

    fmt : str
        A python string format defining the output precision.  Default is .6f,
        which is accurate to 0.0036 arcseconds.

    radunit : str
        This denotes the unit of the radius.

    Returns
    -------
    region_string : str
        ds9 region string

    Examples
    --------
    TODO
    """
    shapelist = to_shape_list(regions, 'DS9', coordsys)
    return shapelist.to_ds9(coordsys, fmt, radunit)


def write_ds9(regions, filename='ds9.reg', coordsys='fk5', fmt='.6f', radunit='deg'):
    """
    Converts list of regions to ds9 string and write to file.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    filename : str
        Filename in which the string is to be written. Default is 'ds9.reg'

    coordsys : str #TODO
        Coordinate system that overrides the coordinate frames of all regions.

    fmt : str
        A python string format defining the output precision.  Default is .6f,
        which is accurate to 0.0036 arcseconds.

    radunit : str
        This denotes the unit of the radius. Default is deg (degrees)
    """
    output = ds9_objects_to_string(regions, coordsys, fmt, radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
