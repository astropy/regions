cimport numpy as np

ctypedef np.float64_t DTYPE_t

cdef int point_in_polygon(double x, double y, np.ndarray[DTYPE_t, ndim=1] vx, np.ndarray[DTYPE_t, ndim=1] vy)
