# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

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

from astropy.tests.pytest_plugins import *

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions
enable_deprecations_as_exceptions()

# The KeyError below is essential in some cases when the package uses other
# astropy affiliated packages.
try:
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    del PYTEST_HEADER_MODULES['h5py']
    del PYTEST_HEADER_MODULES['Pandas']
except KeyError:
    pass

## Uncomment the following lines to display the version number of the
## package rather than the version number of Astropy in the top line when
## running the tests.
import os

## This is to figure out the affiliated package version, rather than
## using Astropy's
try:
    from .version import version
except ImportError:
    version = 'dev'

packagename = os.path.basename(os.path.dirname(__file__))
TESTED_VERSIONS[packagename] = version
