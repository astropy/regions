# Licensed under a 3-clause BSD style license - see LICENSE.rst
from regions.core import PixCoord

parameters = {
    'PixCoord': ['x', 'y'],
}


def pixcoord():
    return PixCoord(1, 2), parameters['PixCoord']
