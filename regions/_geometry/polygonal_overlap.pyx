# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
cimport numpy as np


__all__ = ['polygonal_overlap_grid']


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

from .pnpoly cimport point_in_polygon

__all__ = ['polygonal_overlap_grid']


def polygonal_overlap_grid(double xmin, double xmax, double ymin, double ymax,
                         int nx, int ny,
                         np.ndarray[DTYPE_t, ndim=1] vx,
                         np.ndarray[DTYPE_t, ndim=1] vy,
                         int use_exact,
                         int subpixels):
    """
    polygonal_overlap_grid(xmin, xmax, ymin, ymax, nx, ny, r,
                         use_exact, subpixels)

    Area of overlap between a polygon and a pixel grid.

    Parameters
    ----------
    xmin, xmax, ymin, ymax : float
        Extent of the grid in the x and y direction.
    nx, ny : int
        Grid dimensions.
    vx, vy : `numpy.ndarray`
        The vertices of the polygon
    use_exact : 0 or 1
        If ``1`` calculates exact overlap, if ``0`` uses ``subpixel`` number
        of subpixels to calculate the overlap.
    subpixels : int
        Each pixel resampled by this factor in each dimension, thus each
        pixel is divided into ``subpixels ** 2`` subpixels.

    Returns
    -------
    frac : `~numpy.ndarray` (float)
        2-d array of shape (ny, nx) giving the fraction of the overlap.
    """

    cdef unsigned int i, j
    cdef double x, y, dx, dy, d, pixel_radius
    cdef double bxmin, bxmax, bymin, bymax
    cdef double pxmin, pxcen, pxmax, pymin, pycen, pymax

    # Define output array
    cdef np.ndarray[DTYPE_t, ndim=2] frac = np.zeros([ny, nx], dtype=DTYPE)

    if use_exact == 1:
        raise NotImplementedError("Exact mode has not been implemented for "
                                  "polygonal apertures")

    # Find the width of each element in x and y
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # Define bounding box
    bxmin = vx.min()
    bxmax = vx.max()
    bymin = vy.min()
    bymax = vy.max()

    for i in range(nx):
        pxmin = xmin + i * dx  # lower end of pixel
        pxmax = pxmin + dx  # upper end of pixel
        if pxmax > bxmin and pxmin < bxmax:
            for j in range(ny):
                pymin = ymin + j * dy
                pymax = pymin + dy
                if pymax > bymin and pymin < bymax:
                    frac[j, i] = polygonal_overlap_single_subpixel(pxmin, pymin,
                                                                 pxmax, pymax,
                                                                 vx, vy, subpixels)

    return frac


cdef double polygonal_overlap_single_subpixel(double x0, double y0,
                                            double x1, double y1,
                                            np.ndarray[DTYPE_t, ndim=1] vx,
                                            np.ndarray[DTYPE_t, ndim=1] vy,
                                            int subpixels):
    """
    Return the fraction of overlap between a polygon and a single pixel
    with given extent, using a sub-pixel sampling method.
    """

    cdef unsigned int i, j
    cdef double x, y, dx, dy
    cdef double frac = 0.  # Accumulator.

    dx = (x1 - x0) / subpixels
    dy = (y1 - y0) / subpixels

    x = x0 - 0.5 * dx
    for i in range(subpixels):
        x += dx
        y = y0 - 0.5 * dy
        for j in range(subpixels):
            y += dy
            if point_in_polygon(x, y, vx, vy) == 1:
                frac += 1.

    return frac / (subpixels * subpixels)
