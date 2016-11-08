from distutils.version import LooseVersion

import astropy

ASTROPY_LT_13 = LooseVersion(astropy.__version__) < LooseVersion('1.3')
