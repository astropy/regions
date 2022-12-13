# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    import matplotlib
    HAS_MATPLOTLIB = True
    MPL_VERSION = getattr(matplotlib, '__version__', None)
    if MPL_VERSION is None:
        MPL_VERSION = matplotlib._version.version
    MPL_VERSION = MPL_VERSION.split('.')
    MPL_VERSION = 10 * int(MPL_VERSION[0]) + int(MPL_VERSION[1])
except ImportError:
    HAS_MATPLOTLIB = False
    MPL_VERSION = 0
