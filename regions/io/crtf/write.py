# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from ...core import Region, Regions
from ...core.registry import RegionsRegistry
from ..core import _to_shape_list

__all__ = ['write_crtf', 'crtf_objects_to_string']


@RegionsRegistry.register(Region, 'serialize', 'crtf')
@RegionsRegistry.register(Regions, 'serialize', 'crtf')
def _serialize_crtf(regions, coordsys='fk5', fmt='.6f', radunit='deg'):
    shapelist = _to_shape_list(regions, coordsys)
    return shapelist.to_crtf(coordsys, fmt, radunit)


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
    return _serialize_crtf(regions, coordsys=coordsys, fmt=fmt,
                           radunit=radunit)


@RegionsRegistry.register(Region, 'write', 'crtf')
@RegionsRegistry.register(Regions, 'write', 'crtf')
def write_crtf(regions, filename, coordsys='fk5', fmt='.6f', radunit='deg',
               overwrite=False):
    """
    Convert a list of `~regions.Region` to a CRTF string and write to a
    file.

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

    overwrite : bool, optional
        If True, overwrite the output file if it exists. Raises an
        `OSError` if False and the output file exists. Default is False.
    """
    if os.path.lexists(filename) and not overwrite:
        raise OSError(f'{filename} already exists')

    output = _serialize_crtf(regions, coordsys=coordsys, fmt=fmt,
                             radunit=radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
