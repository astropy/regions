# Licensed under a 3-clause BSD style license - see LICENSE.rst

try:
    import matplotlib
    HAS_MATPLOTLIB = True
    MPL_VER_STR = getattr(matplotlib, '__version__', None)
    if MPL_VER_STR is None:
        MPL_VER_STR = matplotlib._version.version
    MPL_VERSION = MPL_VER_STR.split('.')
    MPL_VERSION = 10 * int(MPL_VERSION[0]) + int(MPL_VERSION[1])
except ImportError:
    HAS_MATPLOTLIB = False
    MPL_VER_STR = 'None'
    MPL_VERSION = 0
