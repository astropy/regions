# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
ASDF converters.
"""
from importlib.util import find_spec


_ASDF_INSTALLED = find_spec('asdf') is not None


__all__ = []


if _ASDF_INSTALLED:
    from .pixel_coords import PixCoordConverter
    __all__ += ['PixCoordConverter',
                ]
