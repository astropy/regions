# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

def trailing_zeros(x):
    """
    Count the number of trailing bit set to zero in ``x``

    Parameters
    ----------
    x : int
        64 bit signed integer

    Returns
    -------
    bits: int
        number of trailing bits set to zero.
    """
    bits = 0
    # convention for x == 0 => return 0
    if x == 0:
        return 0

    if not (x & 0xFFFFFFFF):
        bits += 32
        x >>= 32
    if not (x & 0xFFFF):
        bits += 16
        x >>= 16
    if not (x & 0xFF):
        bits += 8
        x >>= 8
    if not (x & 0xF):
        bits += 4
        x >>= 4
    if not (x & 0x3):
        bits += 2
        x >>= 2
    if not (x & 1):
        bits += 1
        x >>= 1

    return bits
