# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from distutils.version import LooseVersion

try:
    import pytest_arraydiff
except ImportError:
    raise ImportError("The pytest-arraydiff package is required for the tests. "
                      "You can install it with: pip install pytest-arraydiff")
else:
    # We need to remove pytest_arraydiff from the namespace otherwise pytest
    # gets confused, because it tries to interpret pytest_* as a special
    # function name.
    del pytest_arraydiff


from astropy.version import version as astropy_version

if LooseVersion(astropy_version) < LooseVersion('2.0.3'):
    # Astropy is not compatible with the standalone plugins prior this while
    # astroquery requires them, so we need this workaround. This will mess
    # up the test header, but everything else will work.
    from astropy.tests.pytest_plugins import (PYTEST_HEADER_MODULES,
                                              TESTED_VERSIONS)
elif astropy_version < '3.0':
    # With older versions of Astropy, we actually need to import the pytest
    # plugins themselves in order to make them discoverable by pytest.
    from astropy.tests.pytest_plugins import *
else:
    # As of Astropy 3.0, the pytest plugins provided by Astropy are
    # automatically made available when Astropy is installed. This means it's
    # not necessary to import them here, but we still need to import global
    # variables that are used for configuration.
    from astropy.tests.plugins.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS

from astropy.tests.helper import enable_deprecations_as_exceptions


# Uncomment the following line to treat all DeprecationWarnings as
# exceptions
# enable_deprecations_as_exceptions()    # noqa

# Uncomment and customize the following lines to add/remove entries from
# the list of packages for which version numbers are displayed when running
# the tests. Making it pass for KeyError is essential in some cases when
# the package uses other astropy affiliated packages.
try:
    PYTEST_HEADER_MODULES['Cython'] = 'Cython'    # noqa
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'    # noqa
    del PYTEST_HEADER_MODULES['h5py']    # noqa
    del PYTEST_HEADER_MODULES['Pandas']    # noqa
except (NameError, KeyError):  # NameError is needed to support Astropy < 1.0
    pass

# Uncomment the following lines to display the version number of the
# package rather than the version number of Astropy in the top line when
# running the tests.
import os

# This is to figure out the affiliated package version, rather than
# using Astropy's
try:
    from .version import version, astropy_helpers_version
except ImportError:
    version = 'dev'

packagename = os.path.basename(os.path.dirname(__file__))
TESTED_VERSIONS[packagename] = version    # noqa
TESTED_VERSIONS['astropy_helpers'] = astropy_helpers_version
