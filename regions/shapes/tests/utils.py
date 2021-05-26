# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    import matplotlib  # noqa
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
