# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io import registry
from astropy.io.fits.connect import is_fits

from ..core import ShapeList
from .read import read_fits
from .write import write_fits

__all__ = []


registry.register_reader('fits', ShapeList, read_fits)
registry.register_writer('fits', ShapeList, write_fits)
registry.register_identifier('fits', ShapeList, is_fits)
