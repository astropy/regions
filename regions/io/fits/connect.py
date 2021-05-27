# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io import registry
from astropy.io.fits.connect import is_fits

from ..core import ShapeList
from .read import read_fits_region
from .write import write_fits_region

__all__ = []


# FIXME: write_fits_region has its first two arguments backwards as
# compared to write_crtf and write_ds9
def write_fits(regions, filename, *args, **kwargs):
    return write_fits_region(filename, regions, *args, **kwargs)


registry.register_reader('fits', ShapeList, read_fits_region)
registry.register_writer('fits', ShapeList, write_fits)
registry.register_identifier('fits', ShapeList, is_fits)
