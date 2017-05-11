# Licensed under a 3-clause BSD style license - see LICENSE.rst
from distutils.version import LooseVersion

import astropy

ASTROPY_LT_13 = LooseVersion(astropy.__version__) < LooseVersion('1.3')

try:
    import matplotlib  # noqa
    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False

try:
    import shapely  # noqa
    HAS_SHAPELY = True
except:
    HAS_SHAPELY = False
