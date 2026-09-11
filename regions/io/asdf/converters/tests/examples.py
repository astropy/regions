# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from regions.core import PixCoord


parameters = {
    'PixCoord': ['x', 'y'],
}


def pixcoord():
    return PixCoord(1, 2), parameters['PixCoord']


def pixcoord_array():
    return PixCoord(np.array([1, 2]), np.array([3, 4])), parameters['PixCoord']
