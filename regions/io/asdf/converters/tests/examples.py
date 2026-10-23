# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from regions.core import PixCoord

parameters = {
    'PixCoord': ['x', 'y'],
    'CirclePixelRegion': ['center', 'radius', 'meta', 'visual'],
    'CircleSkyRegion': ['center', 'radius', 'meta', 'visual'],
    'CircleSphericalSkyRegion': ['center', 'radius', 'meta', 'visual'],
    'RegionMeta': ['a', 'b'],
}


def pixcoord():
    return PixCoord(1, 2), parameters['PixCoord']


def pixcoord_array():
    return PixCoord(np.array([1, 2]), np.array([3, 4])), parameters['PixCoord']


def circle_pixel_region():
    from regions.shapes.circle import CirclePixelRegion
    return CirclePixelRegion(PixCoord(1, 2), 3), parameters['PixCoord']

def meta():
    from regions.core.metadata import RegionMeta
    return RegionMeta(a=1, b=2), parameters['RegionMeta']
