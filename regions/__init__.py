# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regions is an Astropy coordinated package to provide tools for region
handling.
"""

from ._utils.examples import *  # noqa: F401, F403
from .core import *  # noqa: F401, F403
from .io import *  # noqa: F401, F403
from .shapes import *  # noqa: F401, F403

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
