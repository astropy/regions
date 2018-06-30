# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

from astropy.utils.exceptions import AstropyUserWarning

__all__ = [
    'DS9RegionParserWarning',
    'DS9RegionParserError',
]


class DS9RegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for DS9 region parsing inherited from astropy's
    warnings
    """


class DS9RegionParserError(ValueError):
    """
    A generic error class for DS9 region parsing
    """


# mapping to matplotlib marker symbols, also compatible with CRTF.
valid_symbols_ds9 = dict(circle="o",
                        box="s",
                        diamond="D",
                        x="x",
                        cross="+",
                        arrow="^",
                        boxcircle="*")
