"""
py23_compat.py

Python 2 / 3 compatibility layer.
"""
import sys

__all__ = [
    'int', 'range',
]

PY2 = (sys.version_info.major == 2)

if PY2:
    int = long
    range = xrange
else:
    int = int
    range = range
