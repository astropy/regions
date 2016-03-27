# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Functions for performing aperture photometry on 2-D arrays."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np


def get_phot_extents(shape, positions, extents):
    """
    Get the photometry extents and check if the apertures is fully out of data.

    Parameters
    ----------
    shape : tuple
        The shape tuple describing an array

    Returns
    -------
    extents : dict
        The ``extents`` dictionary contains 3 elements:

        * ``'ood_filter'``
            A boolean array with `True` elements where the aperture is
            falling out of the data region.
        * ``'pixel_extent'``
            x_min, x_max, y_min, y_max : Refined extent of apertures with
            data coverage.
        * ``'phot_extent'``
            x_pmin, x_pmax, y_pmin, y_pmax: Extent centered to the 0, 0
            positions as required by the `~photutils.geometry` functions.
    """

    # Check if an aperture is fully out of data
    ood_filter = np.logical_or(extents[:, 0] >= shape[1],
                               extents[:, 1] <= 0)
    np.logical_or(ood_filter, extents[:, 2] >= shape[0],
                  out=ood_filter)
    np.logical_or(ood_filter, extents[:, 3] <= 0, out=ood_filter)

    # TODO check whether it makes sense to have negative pixel
    # coordinate, one could imagine a stackes image where the reference
    # was a bit offset from some of the images? Or in those cases just
    # give Skycoord to the Aperture and it should deal with the
    # conversion for the actual case?
    x_min = np.maximum(extents[:, 0], 0)
    x_max = np.minimum(extents[:, 1], shape[1])
    y_min = np.maximum(extents[:, 2], 0)
    y_max = np.minimum(extents[:, 3], shape[0])

    x_pmin = x_min - positions[:, 0] - 0.5
    x_pmax = x_max - positions[:, 0] - 0.5
    y_pmin = y_min - positions[:, 1] - 0.5
    y_pmax = y_max - positions[:, 1] - 0.5

    # TODO: check whether any pixel is nan in data[y_min[i]:y_max[i],
    # x_min[i]:x_max[i])), if yes return something valid rather than nan

    pixel_extent = [x_min, x_max, y_min, y_max]
    phot_extent = [x_pmin, x_pmax, y_pmin, y_pmax]

    return ood_filter, pixel_extent, phot_extent
