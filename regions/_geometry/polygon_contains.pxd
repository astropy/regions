# Licensed under a 3-clause BSD style license - see LICENSE.rst
# cython: language_level=3
"""
Declarations needed to cimport the point-in-polygon functions into other
Cython files. All functions are pure C math functions that are safe to
call without the GIL.
"""

cdef int point_in_polygon(double x, double y, const double *vx,
                          const double *vy, int n) noexcept nogil
cdef int point_in_polygon_covers(double x, double y, const double *vx,
                                 const double *vy, int n) noexcept nogil
cdef int point_in_polygon_even_odd(double x, double y, const double *vx,
                                   const double *vy, int n) noexcept nogil
