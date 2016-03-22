# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an experimental package for representing regions
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    pass

from .io.read_ds9 import parse_ds9
