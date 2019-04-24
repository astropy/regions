# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

from astropy.utils.exceptions import AstropyUserWarning

__all__ = [
    'FITSRegionParserWarning',
    'FITSRegionParserError',
]


class FITSRegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for FITS region parsing
    """


class FITSRegionParserError(ValueError):
    """
    A generic error class for FITS region parsing
    """


language_spec = {'CIRCLE': ['X0', 'Y0', 'R0'],
                 'POINT': ['X0', 'Y0'],
                 'BOX': ['X0', 'Y0', 'R0', 'R1'],
                 'ANNULUS': ['X0', 'Y0', 'R0', 'R1'],
                 'ELLIPSE': ['X0', 'Y0', 'R0', 'R1', 'ROTANG0'],
                 'ELLIPTANNULUS': ['X0', 'Y0', 'R0', 'R1', 'R2', 'R3', 'ROTANG0'],
                 'ROTBOX': ['X0', 'Y0', 'R0', 'R1', 'ROTANG0'],
                 'RECTANGLE': ['X0', 'X1', 'Y0', 'Y1'],
                 'ROTRECTANGLE': ['X0', 'X1', 'Y0', 'Y1', 'ROTANG0'],
                 'POLYGON': ['X', 'Y'],
                 'PIE': ['X0', 'Y0', 'ROTANG0', 'ROTANG1']
                 }
language_spec['SECTOR'] = language_spec['PIE']
