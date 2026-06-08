# Licensed under a 3-clause BSD style license - see LICENSE.rst
# cython: language_level=3
#
# The code in this file was adapted from code written by Greg von Winckel:
#
# https://github.com/gregvw/pnpoly
#
# and which was released under the following license:
#
# ----------------------------------------------------------------------------
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Greg von Winckel
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ----------------------------------------------------------------------------
#
# This code was itself adapted from code written by W. Randolph Franklin:
#
# http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
#
# and released under the following license:
#
# Copyright (c) 1970-2003, Wm. Randolph Franklin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimers.
#
# Redistributions in binary form must reproduce the above copyright notice in
# the documentation and/or other materials provided with the distribution.
#
# The name of W. Randolph Franklin may not be used to endorse or promote
# products derived from this Software without specific prior written permission.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
cimport numpy as np

DTYPE_BOOL = bool
ctypedef np.uint8_t DTYPE_BOOL_t


def points_in_polygon(np.ndarray[DTYPE_t, ndim=1] x,
                      np.ndarray[DTYPE_t, ndim=1] y,
                      np.ndarray[DTYPE_t, ndim=1] vx,
                      np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether points are strictly inside a polygon (excluding
    boundary).

    This is consistent with Shapely's ``contains`` function and DE-9IM
    semantics.

    Parameters
    ----------
    x, y : ndarray
        The x and y coordinates of the test points.
    vx, vy : ndarray
        The x and y coordinates of the polygon vertices.

    Returns
    -------
    result : ndarray of bool
        True if the point is strictly inside the polygon, False otherwise.
    """
    cdef int i, n
    cdef np.ndarray[np.uint8_t, ndim=1, cast=True] result

    n = x.shape[0]

    result = np.zeros(n, DTYPE_BOOL)

    for i in range(n):
        result[i] = point_in_polygon(x[i], y[i], vx, vy)

    return result


def points_in_polygon_covers(np.ndarray[DTYPE_t, ndim=1] x,
                             np.ndarray[DTYPE_t, ndim=1] y,
                             np.ndarray[DTYPE_t, ndim=1] vx,
                             np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether points are inside or on the boundary of a polygon.

    This is consistent with Shapely's ``covers`` function and DE-9IM
    semantics.

    Parameters
    ----------
    x, y : ndarray
        The x and y coordinates of the test points.
    vx, vy : ndarray
        The x and y coordinates of the polygon vertices.

    Returns
    -------
    result : ndarray of bool
        True if the point is inside or on the boundary of the polygon,
        False otherwise.
    """
    cdef int i, n
    cdef np.ndarray[np.uint8_t, ndim=1, cast=True] result

    n = x.shape[0]

    result = np.zeros(n, DTYPE_BOOL)

    for i in range(n):
        result[i] = point_in_polygon_covers(x[i], y[i], vx, vy)

    return result


cdef inline int point_on_segment(double px, double py,
                                  double x1, double y1,
                                  double x2, double y2) noexcept nogil:
    """
    Check whether the point (px, py) lies on the line segment from
    (x1, y1) to (x2, y2).

    The collinearity tolerance is expressed as a perpendicular
    distance from the point to the segment so that the test behaves
    consistently regardless of the segment length or the magnitude of
    the coordinates.

    The tolerance is an absolute distance in pixel units (1e-10 px), not
    a relative tolerance, so it is most appropriate for the pixel-scale
    coordinates used by the region classes.
    """
    cdef double x_min, x_max, y_min, y_max, cross, seg_len_sq
    cdef double tol = 1e-10  # perpendicular distance tolerance (pixels)

    # Check bounding box
    if x1 < x2:
        x_min = x1
        x_max = x2
    else:
        x_min = x2
        x_max = x1

    if y1 < y2:
        y_min = y1
        y_max = y2
    else:
        y_min = y2
        y_max = y1

    if px < x_min or px > x_max or py < y_min or py > y_max:
        return 0

    seg_len_sq = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)
    if seg_len_sq == 0:
        # Degenerate segment (duplicate vertices); treat it as a point.
        return 1 if (px == x1 and py == y1) else 0

    # The cross product magnitude equals perpendicular_distance times
    # segment_length, so comparing its square to tol**2 * seg_len_sq
    # tests the perpendicular distance against tol independently of the
    # segment length.
    cross = (px - x1) * (y2 - y1) - (py - y1) * (x2 - x1)
    if cross * cross > tol * tol * seg_len_sq:
        return 0

    return 1


cdef int point_in_polygon_even_odd(double x, double y,
                                   np.ndarray[DTYPE_t, ndim=1] vx,
                                   np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether the point (x, y) is inside a polygon using the
    even-odd (ray-casting) rule, ignoring the polygon boundary.

    Boundary behavior is undefined here; callers that need well-defined
    boundary semantics should use `point_in_polygon` (boundary excluded)
    or `point_in_polygon_covers` (boundary included).

    This uses the even-odd rule, as described here:

        https://en.wikipedia.org/wiki/Even–odd_rule
    """
    cdef int i, j, n
    cdef int result = 0

    n = vx.shape[0]

    for i in range(n):
        j = (i + n - 1) % n
        if(((vy[i] > y) != (vy[j] > y)) and
            (x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i])):
            result += 1

    return result % 2


cdef int point_in_polygon(double x, double y,
                          np.ndarray[DTYPE_t, ndim=1] vx,
                          np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether a test point (x, y) is strictly inside a polygon
    defined by a set of vertices (vx, vy).

    Points on the boundary (edges and vertices) are NOT considered
    inside, consistent with Shapely's ``contains`` function and DE-9IM
    semantics.
    """
    cdef int i, j, n

    # The even-odd (ray-casting) rule is reliable for points that are
    # strictly inside or strictly outside; it is only ambiguous on the
    # boundary. Points reported as outside are not inside and are not
    # strictly inside even if they lie on an edge, so no edge scan is
    # needed for them. Only points reported as inside require an edge
    # scan to exclude boundary points.
    if point_in_polygon_even_odd(x, y, vx, vy) == 0:
        return 0

    n = vx.shape[0]
    for i in range(n):
        j = (i + n - 1) % n
        if point_on_segment(x, y, vx[i], vy[i], vx[j], vy[j]):
            return 0

    return 1


cdef int point_in_polygon_covers(double x, double y,
                                 np.ndarray[DTYPE_t, ndim=1] vx,
                                 np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether a test point (x, y) is inside or on the boundary
    of a polygon defined by a set of vertices (vx, vy).

    Points on the boundary (edges and vertices) ARE considered inside,
    consistent with Shapely's ``covers`` function and DE-9IM semantics.
    """
    cdef int i, j, n

    # The even-odd (ray-casting) rule is reliable for points that are
    # strictly inside or strictly outside; it is only ambiguous on the
    # boundary. Points reported as inside are covered regardless of the
    # boundary, so no edge scan is needed for them. Only points reported
    # as outside require an edge scan to include boundary points.
    if point_in_polygon_even_odd(x, y, vx, vy) == 1:
        return 1

    n = vx.shape[0]
    for i in range(n):
        j = (i + n - 1) % n
        if point_on_segment(x, y, vx[i], vy[i], vx[j], vy[j]):
            return 1

    return 0
