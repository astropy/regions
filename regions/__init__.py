# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an in-development package for region handling based on Astropy.

The goal is to merge the functionality from pyregion and photutils apertures
and then after some time propose this package for inclusion in the Astropy core.

* Code : https://github.com/astropy/regions
* Docs : http://astropy-regions.readthedocs.io/en/latest/
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from ._utils.examples import *
    from .core import *
    from .shapes import *
    from .io import *
