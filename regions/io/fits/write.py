# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.io import fits

from ..core import to_shape_list, SkyRegion

__all__ = ['write_fits_region', 'fits_region_objects_to_table']


def fits_region_objects_to_table(regions):
    """
    Convert a list of `~regions.Region` objects to a FITS region table.

    See :ref:`gs-ds9`

    Parameters
    ----------
    regions : list
        A list of `regions.Region` objects.

    Returns
    -------
    region_string : `~astropy.table.Table`
        A FITS region table.
    """
    for reg in regions:
        if isinstance(reg, SkyRegion):
            raise TypeError('Every region must be a pixel region')

    shape_list = to_shape_list(regions, coordinate_system='image')
    return shape_list.to_fits()


def write_fits_region(filename, regions, header=None, overwrite=False):
    """
    Convert a list of `~regions.Region` to a FITS region table and write
    to a file.

    See :ref:`gs-fits`

    Parameters
    ----------
    filename : str
        The filename in which the table is to be written.

    regions : list
        A list of `~regions.Region` objects.

    header : `~astropy.io.fits.Header`, optional
        The FITS header.

    overwrite : bool, optional
        If True, overwrite the output file if it exists. Raises an
        `OSError` if False and the output file exists. Default is False.
    """
    output = fits_region_objects_to_table(regions)
    bin_table = fits.BinTableHDU(data=output, header=header)
    bin_table.writeto(filename, overwrite=overwrite)
