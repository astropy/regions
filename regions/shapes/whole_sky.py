# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines a whole-sky spherical sky region, for use with
compound region logic and over-the-pole logic.
"""

import numpy as np

from regions.core.core import SphericalSkyRegion


class WholeSphericalSkyRegion(SphericalSkyRegion):
    """
    Spherical region representing the whole sky.

    Implemented to handle compound region logic, particularly with
    ranges and over-the-pole pole logic.
    """

    def contains(self, coord):
        if coord.isscalar:
            return True
        else:
            return np.ones(coord.shape, dtype=bool)

    @property
    def bounding_circle(self):
        # Not defined
        return None

    @property
    def bounding_lonlat(self):
        # Not defined
        return None

    def transform_to(self, frame):
        # TODO: handle offset origin transformations

        return self.__class__()

    def to_sky(
        self,
        wcs=None,
        include_boundary_distortions=False,
        discretize_kwargs=None
    ):
        # Not defined
        raise NotImplementedError('Not defined.')

    def to_pixel(
        self,
        wcs=None,
        include_boundary_distortions=False,
        discretize_kwargs=None
    ):
        # Not defined
        raise NotImplementedError('Not defined.')
