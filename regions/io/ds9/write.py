# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy import coordinates
from ... import shapes
from ..core import to_shape_list

__all__ = [
    'write_ds9',
    'ds9_objects_to_string',
]

coordsys_name_mapping = dict(zip(coordinates.frame_transform_graph.get_names(),
                                 coordinates.frame_transform_graph.get_names()))
coordsys_name_mapping['ecliptic'] = 'geocentrictrueecliptic'  # needs expert attention TODO


def ds9_objects_to_string(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    """Convert list of regions to ds9 region strings.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    Returns
    -------
    region_string : str
        ds9 region string
    fmt : str
        A python string format defining the output precision.  Default is .6f,
        which is accurate to 0.0036 arcseconds.


    Examples
    --------
    TODO
    """
    shapelist = to_shape_list(regions, 'DS9', coordsys)
    return shapelist.to_ds9(fmt, radunit, coordsys)


def write_ds9(regions, filename='ds9.reg', coordsys='fk5'):
    """
    Convert list of regions to ds9 string and write to file.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects
    filename : str
        Filename
    coordsys : {TODO}
        Coordinate system
    """
    output = ds9_objects_to_string(regions, coordsys)
    with open(filename, 'w') as fh:
        fh.write(output)
