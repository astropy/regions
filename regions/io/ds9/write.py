# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from ...core import Region, Regions
from ...core.registry import RegionsRegistry
from ..core import _to_shape_list

__all__ = []


@RegionsRegistry.register(Region, 'serialize', 'ds9')
@RegionsRegistry.register(Regions, 'serialize', 'ds9')
def _serialize_ds9(regions, fmt='.6f', radunit='deg'):
    shapelist = _to_shape_list(regions)
    return shapelist.to_ds9(fmt, radunit)


@RegionsRegistry.register(Region, 'write', 'ds9')
@RegionsRegistry.register(Regions, 'write', 'ds9')
def _write_ds9(regions, filename, fmt='.6f', radunit='deg',
               overwrite=False):
    """
    Convert a list of `~regions.Region` to a DS9 string and write to a
    file.

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    filename : str
        The filename in which the string is to be written.

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

    output = _serialize_ds9(regions, fmt=fmt, radunit=radunit)
    with open(filename, 'w') as fh:
        fh.write(output)
