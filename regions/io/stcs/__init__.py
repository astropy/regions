# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage provides tools for reading and writing STC-S region files.

The Space-Time Coordinate (STC) string representation (STC-S) is an IVOA
standard for describing spatial and temporal regions and coordinates.
"""

from .read import *  # noqa: F401, F403
from .write import *  # noqa: F401, F403
from .core import *  # noqa: F401, F403
from .connect import *  # noqa: F401, F403
