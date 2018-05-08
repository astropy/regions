# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

from astropy.utils.exceptions import AstropyUserWarning

__all__ = [
    'CRTFRegionParserWarning',
    'CRTFRegionParserError',
]


class CRTFRegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for CRTF region parsing
    """


class CRTFRegionParserError(ValueError):
    """
    A generic error class for CRTF region parsing
    """
