# Licensed under a 3-clause BSD style license - see LICENSE.rst
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

DTYPE_BOOL = np.bool
ctypedef np.uint8_t DTYPE_BOOL_t


def points_in_polygon(np.ndarray[DTYPE_t, ndim=1] x,
                      np.ndarray[DTYPE_t, ndim=1] y,
                      np.ndarray[DTYPE_t, ndim=1] vx,
                      np.ndarray[DTYPE_t, ndim=1] vy):

    cdef int i, n
    cdef np.ndarray[np.uint8_t, ndim=1, cast=True] result

    n = x.shape[0]

    result = np.zeros(n, DTYPE_BOOL)

    for i in range(n):
        result[i] = point_in_polygon(x[i], y[i], vx, vy)

    return result


cdef int point_in_polygon(double x, double y,
                          np.ndarray[DTYPE_t, ndim=1] vx,
                          np.ndarray[DTYPE_t, ndim=1] vy):
    """
    Determine whether a test point (x, y) is within a polygon defined by a set
    of vertices (vx, vy).

    This uses the even-odd rule, as described here:

        https://en.wikipedia.org/wiki/Even–odd_rule
    """

    cdef int i, j, k, m, n
    cdef int result

    n = vx.shape[0]

    result = 0

    for i in range(n):
        j = (i + n - 1) % n
        if(((vy[i] > y) != (vy[j] > y)) and
            (x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i])):
            result += 1

    return result % 2
