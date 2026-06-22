# Licensed under a 3-clause BSD style license - see LICENSE.rst
# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# cython: freethreading_compatible=True
"""
Tools for testing whether points are inside polygons.

The public API is provided by the Python functions `points_in_polygon`
and `points_in_polygon_covers`, which operate on arrays of points and
return boolean arrays indicating whether each point is inside the
polygon (excluding or including the boundary, respectively).

The cdef functions are not intended to be called from Python code. They
are pure C math functions declared ``noexcept nogil`` so they can be
called without the GIL, including from multiple threads on free-threaded
Python builds. Their signatures are exported via polygon_contains.pxd.
"""

import numpy as np

cimport numpy as np

__all__ = ['points_in_polygon', 'points_in_polygon_covers']


DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


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
    cdef Py_ssize_t i, n
    cdef int n_poly
    cdef const double[::1] x_view = np.ascontiguousarray(x)
    cdef const double[::1] y_view = np.ascontiguousarray(y)
    cdef const double[::1] vx_view = np.ascontiguousarray(vx)
    cdef const double[::1] vy_view = np.ascontiguousarray(vy)
    cdef const double *vx_ptr
    cdef const double *vy_ptr

    n = x_view.shape[0]
    n_poly = vx_view.shape[0]

    cdef np.ndarray[np.uint8_t, ndim=1] result = np.zeros(n, dtype=np.uint8)
    cdef np.uint8_t[::1] result_view = result

    if n_poly > 0:
        vx_ptr = &vx_view[0]
        vy_ptr = &vy_view[0]
        with nogil:
            for i in range(n):
                result_view[i] = point_in_polygon(x_view[i], y_view[i],
                                                  vx_ptr, vy_ptr, n_poly)

    return result.astype(bool)


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
    cdef Py_ssize_t i, n
    cdef int n_poly
    cdef const double[::1] x_view = np.ascontiguousarray(x)
    cdef const double[::1] y_view = np.ascontiguousarray(y)
    cdef const double[::1] vx_view = np.ascontiguousarray(vx)
    cdef const double[::1] vy_view = np.ascontiguousarray(vy)
    cdef const double *vx_ptr
    cdef const double *vy_ptr

    n = x_view.shape[0]
    n_poly = vx_view.shape[0]

    cdef np.ndarray[np.uint8_t, ndim=1] result = np.zeros(n, dtype=np.uint8)
    cdef np.uint8_t[::1] result_view = result

    if n_poly > 0:
        vx_ptr = &vx_view[0]
        vy_ptr = &vy_view[0]
        with nogil:
            for i in range(n):
                result_view[i] = point_in_polygon_covers(x_view[i], y_view[i],
                                                         vx_ptr, vy_ptr, n_poly)

    return result.astype(bool)


cdef int point_in_polygon(double x, double y, const double *vx,
                          const double *vy, int n) noexcept nogil:
    """
    Determine whether a test point (x, y) is strictly inside a polygon
    defined by a set of vertices (vx, vy).

    Points on the boundary (edges and vertices) are NOT considered
    inside, consistent with Shapely's ``contains`` function and DE-9IM
    semantics.
    """
    cdef int i, j

    # The even-odd (ray-casting) rule is reliable for points that are
    # strictly inside or strictly outside; it is only ambiguous on the
    # boundary. Points reported as outside are not inside and are not
    # strictly inside even if they lie on an edge, so no edge scan is
    # needed for them. Only points reported as inside require an edge
    # scan to exclude boundary points.
    if point_in_polygon_even_odd(x, y, vx, vy, n) == 0:
        return 0

    for i in range(n):
        j = (i + n - 1) % n
        if point_on_segment(x, y, vx[i], vy[i], vx[j], vy[j]):
            return 0

    return 1


cdef int point_in_polygon_covers(double x, double y, const double *vx,
                                 const double *vy, int n) noexcept nogil:
    """
    Determine whether a test point (x, y) is inside or on the boundary
    of a polygon defined by a set of vertices (vx, vy).

    Points on the boundary (edges and vertices) ARE considered inside,
    consistent with Shapely's ``covers`` function and DE-9IM semantics.
    """
    cdef int i, j

    # The even-odd (ray-casting) rule is reliable for points that are
    # strictly inside or strictly outside; it is only ambiguous on the
    # boundary. Points reported as inside are covered regardless of the
    # boundary, so no edge scan is needed for them. Only points reported
    # as outside require an edge scan to include boundary points.
    if point_in_polygon_even_odd(x, y, vx, vy, n) == 1:
        return 1

    for i in range(n):
        j = (i + n - 1) % n
        if point_on_segment(x, y, vx[i], vy[i], vx[j], vy[j]):
            return 1

    return 0


cdef int point_in_polygon_even_odd(double x, double y, const double *vx,
                                   const double *vy, int n) noexcept nogil:
    """
    Determine whether the point (x, y) is inside a polygon using the
    even-odd (ray-casting) rule, ignoring the polygon boundary.

    Boundary behavior is undefined here; callers that need well-defined
    boundary semantics should use `point_in_polygon` (boundary excluded)
    or `point_in_polygon_covers` (boundary included).

    This uses the even-odd rule, as described here:

        https://en.wikipedia.org/wiki/Even-odd_rule
    """
    cdef int i, j
    cdef int result = 0

    for i in range(n):
        j = (i + n - 1) % n
        if (((vy[i] > y) != (vy[j] > y)) and
                (x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i])):
            result += 1

    return result % 2


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
