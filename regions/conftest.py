# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure. It needs to live inside the package in order for it to
# get picked up when running the tests inside an interpreter using
# packagename.test

try:
    import pytest_arraydiff
except ImportError:
    raise ImportError("The pytest-arraydiff package is required to run the "
                      "tests. You can install it with: pip install "
                      "pytest-arraydiff.")
else:
    # We need to remove pytest_arraydiff from the namespace otherwise pytest
    # gets confused, because it tries to interpret pytest_* as a special
    # function name.
    del pytest_arraydiff

try:
    from pytest_astropy_header.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):
    if ASTROPY_HEADER:
        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the
        # list of packages for which version numbers are displayed when
        # running the tests.
        PYTEST_HEADER_MODULES['Cython'] = 'Cython'
        PYTEST_HEADER_MODULES['Numpy'] = 'numpy'
        PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
        PYTEST_HEADER_MODULES['Matplotlib'] = 'matplotlib'
        PYTEST_HEADER_MODULES['Shapely'] = 'shapely'
        PYTEST_HEADER_MODULES.pop('scipy', None)
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES.pop('h5py', None)

        from regions import __version__
        TESTED_VERSIONS['regions'] = __version__
