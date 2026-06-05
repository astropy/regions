# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Configuration file for the pytest test suite.
"""

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
        PYTEST_HEADER_MODULES.clear()
        deps = ['NumPy', 'Matplotlib', 'Astropy', 'Shapely']
        for dep in deps:
            PYTEST_HEADER_MODULES[dep] = dep.lower()

        from regions import __version__
        TESTED_VERSIONS['regions'] = __version__
