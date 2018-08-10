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

    >>> from regions import CirclePixelRegion, PixCoord
    >>> reg_pixel = CirclePixelRegion(PixCoord(1, 2), 5)
    >>> table = fits_region_objects_to_table([reg_pixel])
    >>> print(table)
    X [1] Y [1] SHAPE    R [4]    ROTANG COMPONENT
    pix   pix            pix      deg
    ----- ----- ------ ---------- ------ ---------
      1.0   2.0 circle 5.0 .. 0.0      0         1

    """
    for reg in regions:
        if isinstance(reg, SkyRegion):
            raise TypeError('Every region must be a pixel region'.format(reg))

    shape_list = to_shape_list(regions, coordinate_system='image')
    return shape_list.to_fits()


def write_fits_region(filename, regions, header=None):
    """
    Converts list of regions to FITS region table and write to a file.

    Parameters
    ----------
    filename: str
        Filename in which the table is to be written. Default is 'new.fits'
    regions: list
        List of `regions.Region` objects
    header: `~astropy.io.fits.header.Header` object
        The FITS header.

    Examples
    --------

    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from astropy.io import fits
    >>> file_sample = get_pkg_data_filename('data/region.fits', package='regions.io.fits.tests')
    >>> from regions import CirclePixelRegion, PixCoord, write_fits_region
    >>> reg_pixel = CirclePixelRegion(PixCoord(1, 2), 5)
    >>> hdul = fits.open(filename)
    >>> write_fits_region('region_output.fits', regions=[reg_pixel], header=hdul[1].header)

    """
    output = fits_region_objects_to_table(regions)

    bin_table = fits.BinTableHDU(data=output, header=header)
    bin_table.writeto(filename)
