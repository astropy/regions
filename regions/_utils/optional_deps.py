# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Checks for optional dependencies using lazy import from `PEP 562
<https://www.python.org/dev/peps/pep-0562/>`_.
"""

import importlib

# This list is a duplicate of the dependencies in pyproject.toml "all"
# and "test". gwcs is used only in tests.
optional_deps = ['gwcs', 'matplotlib', 'shapely']
deps = {key.upper(): key for key in optional_deps}
__all__ = [f'HAS_{pkg}' for pkg in deps]


def __getattr__(name):
    if name in __all__:
        try:
            importlib.import_module(deps[name[4:]])
        except ImportError:
            return False
        return True

    raise AttributeError(f'Module {__name__!r} has no attribute {name!r}.')
