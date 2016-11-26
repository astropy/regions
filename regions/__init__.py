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
    from .utils.examples import *
    from .core import *
    from .io import *
    from .shapes import *
else:
    import sys
    if 'test' in sys.argv:
        try:
            import pytest_arraydiff
        except ImportError:
            raise ImportError("The pytest-arraydiff package is required for the tests. "
                              "You can install it with: pip install pytest-arraydiff")
        else:
            # Just make sure the plugin import works to avoid obscure errors later
            import pytest_arraydiff.plugin
