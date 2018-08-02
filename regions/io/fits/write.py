# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from astropy.io import fits

from ..core import to_shape_list, SkyRegion


__all__ = [
    'write_fits_region',
    'fits_region_objects_to_table',
]


def fits_region_objects_to_table(regions):
    """
    Converts list of regions to FITS region table.

    Parameters
    ----------
    regions : list
        List of `regions.Region` objects

    Returns
    -------
    region_string : `~astropy.table.Table`
        FITS region table

    Examples
    --------
    TODO
    """
    for reg in regions:
        if isinstance(reg, SkyRegion):
            raise TypeError('This {} must be a pixel region'.format(reg))

    shape_list = to_shape_list(regions, coordinate_system='image')
    return shape_list.to_fits()


def write_fits_region(regions, header=None, filename='new.fits'):
    """
    Converts list of regions to FITS region table and write to a file.

    Parameters
    ----------
    regions: list
        List of `regions.Region` objects

    header: `~astropy.io.fits.header.Header` object
        The FITS header.

    filename: str
        Filename in which the table is to be written. Default is 'new.fits'
    """
    output = fits_region_objects_to_table(regions)

    bin_table = fits.BinTableHDU(data=output, header=header)
    bin_table.writeto(filename)
