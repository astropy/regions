# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from .py23_compat import range, int

import copy
import numpy as np

from astropy_healpix import uniq_to_level_ipix

class IntervalSet:
    """
    Internal data structure for storing a MOC.

    Please refer to the doc example page. See :ref:`moc`.

    Parameters
    ----------
    data : `~numpy.ndarray`
        A N rows by 2 columns numpy array, each row representing one interval.
    make_consistent : bool, optional
        True by default. Remove the overlapping intervals from the data.
    """
    HPY_MAX_DEPTH = 29

    def __init__(self, data=None, make_consistent=True):
        data = np.array([], dtype=np.int64) if data is None else data
        self._data = data
        # Merging of the overlapping intervals
        if make_consistent:
            self._merge_intervals()

    def copy(self):
        """
        Deepcopy of self.

        Returns
        -------
        itv_s : `IntervalSet`
            A copy of self
        """
        return copy.deepcopy(self)

    def __repr__(self):
        # Simply print the numpy array storing the intervals
        return "{0}".format(self._data)

    def __eq__(self, other):
        """
        Test equality of two `IntervalSet` objects

        Parameters
        ----------
        other : `IntervalSet`
            The other `IntervalSet` object

        Returns
        -------
        res : bool
            True if all the intervals from the two `IntervalSet` are equal one by one.
        """
        return np.all(self._data == other._data)

    @property
    def min(self):
        return self._data.min()

    @property
    def max(self):
        return self._data.max()

    def empty(self):
        """
        Return True if there is no interval stored.
        """
        return self._data.size == 0

    def _merge_intervals(self):
        """
        Merge overlapping intervals.

        This method is called only one time at the end of the constructor.
        """
        res = []
        start = stop = None
        itvs_l = self._data.tolist()
        
        for itv in sorted(itvs_l):
            if start is None:
                start, stop = itv
                continue

            # gap between intervals
            if itv[0] > stop:
                res.append((start, stop))
                start, stop = itv
            else:
                # merge intervals
                if itv[1] > stop:
                    stop = itv[1]

        if start is not None and stop is not None:
            res.append((start, stop))

        self._data = np.asarray(res, dtype=np.int64)

    def union(self, other):
        """
        Returns the union between two `IntervalSet` instances.

        Parameters
        ----------
        other : `IntervalSet`
            The other `IntervalSet` instance.
        Returns
        -------
        res : `IntervalSet`
            A new `IntervalSet` instance.
        """
        res = IntervalSet()
        if other.empty():
            res._data = self._data
        elif self.empty():
            res._data = other._data
        else:
            # res has no overlapping intervals
            res._data = IntervalSet.merge(self._data,
                                          other._data,
                                          lambda in_a, in_b: in_a or in_b)
        return res

    def difference(self, other):
        """
        Returns the difference between two `IntervalSet` instances.

        Parameters
        ----------
        other : `IntervalSet`
            The other `IntervalSet` instance.
        Returns
        -------
        res : `IntervalSet`
            A new `IntervalSet` instance.
        """
        res = IntervalSet()
        if other.empty() or self.empty():
            res._data = self._data
        else:
            res._data = IntervalSet.merge(self._data,
                                          other._data,
                                          lambda in_a, in_b: in_a and not in_b)
        return res

    def intersection(self, other):
        """
        Return the intersection between self and ``another_is``.

        Parameters
        ----------
        other : `IntervalSet`
            The other `IntervalSet` object.
        Returns
        -------
        res : `IntervalSet`
            A new `IntervalSet` instance.
        """
        res = IntervalSet()
        if not other.empty() and not self.empty():
            res._data = IntervalSet.merge(self._data,
                                          other._data,
                                          lambda in_a, in_b: in_a and in_b)
        return res

    @classmethod
    def to_uniq_itv_s(cls, itv_s):
        """
        Convert invervals (i.e. using the nested numbering scheme) to intervals of uniq numbers.

        Parameters
        ----------
        itv_s : `IntervalSet`
            A `IntervalSet` instance.

        Returns
        -------
        result : `IntervalSet`
            A new `IntervalSet` instance containing intervals in the uniq format.
        """
        r2 = itv_s.copy()
        uniq_itvs_l = []

        if r2.empty():
            return IntervalSet()

        depth = 0
        while not r2.empty():
            shift = int(2 * (IntervalSet.HPY_MAX_DEPTH - depth))
            ofs = (int(1) << shift) - 1
            ofs2 = int(1) << (2 * depth + 2)

            r4 = []
            for iv in r2._data:
                a = (int(iv[0]) + ofs) >> shift
                b = int(iv[1]) >> shift

                c = a << shift
                d = b << shift
                if d > c:
                    r4.append((c, d))
                    uniq_itvs_l.append((a + ofs2, b + ofs2))

            if len(r4) > 0:
                r4_is = IntervalSet(np.asarray(r4, dtype=np.int64))
                r2 = r2.difference(r4_is)

            depth += 1

        uniq_itvs = np.asarray(uniq_itvs_l, dtype=np.int64)
        return IntervalSet(uniq_itvs)

    @classmethod
    def from_uniq_itv_s(cls, uniq_itv_s):
        """
        Convert uniq numbers intervals to nested intervals.

        Parameters
        ----------
        uniq_itv_s : `IntervalSet`
            A uniq formatted `IntervalSet` instance.

        Returns
        -------
        itv_s : `IntervalSet`
            A `IntervalSet` instance containing nested indices.
        """
        itv_s = IntervalSet()
        # Appending a list is faster than appending a numpy array
        rtmp = []
        last_level = 0
        uniq_itvs = uniq_itv_s._data
        shift = 2 * IntervalSet.HPY_MAX_DEPTH
        for uniq_itv in uniq_itvs:
            for j in np.arange(uniq_itv[0], uniq_itv[1], dtype=np.int64):
                level, ipix = uniq_to_level_ipix(j)

                if level != last_level:
                    itv_s = itv_s.union(IntervalSet(np.asarray(rtmp, dtype=np.int64)))
                    rtmp = []
                    last_level = level
                    shift = 2 * (IntervalSet.HPY_MAX_DEPTH - level)

                rtmp.append((ipix << shift, (ipix + 1) << shift))

        itv_s = itv_s.union(IntervalSet(np.asarray(rtmp, dtype=np.int64)))
        return itv_s

    @staticmethod
    def merge(a_itvs, b_itvs, op):
        """
        Merge two lists of intervals according to the boolean function op

        ``a_itvs`` and ``b_itvs`` need to be sorted and consistent (i.e. no overlapping intervals).
        This operation keeps the resulting `IntervalSet` sorted and consistent.

        Parameters
        ----------
        a_itvs : `~numpy.ndarray`
            A N rows by 2 columns numpy array storing non overlapping and sorted intervals.
        b_itvs : `~numpy.ndarray`
            A N rows by 2 columns numpy array storing non overlapping and sorted intervals.
        op : `function`
            Lambda function taking two params and returning the result of the operation between
            these two params.
            Exemple : lambda a, b: a and b describes the intersection of ``a_itvs`` and
            ``b_itvs``.

        Returns
        -------
        res : `~numpy.ndarray`
            A N rows by 2 columns numpy array containing intervals resulting from the op between ``a_itvs`` and ``b_itvs``.
        """
        a_endpoints = a_itvs.flatten().tolist()
        b_endpoints = b_itvs.flatten().tolist()

        sentinel = max(a_endpoints[-1], b_endpoints[-1]) + 1

        a_endpoints += [sentinel]
        b_endpoints += [sentinel]

        a_index = 0
        b_index = 0

        res = []

        scan = min(a_endpoints[0], b_endpoints[0])
        while scan < sentinel:
            in_a = not ((scan < a_endpoints[a_index]) ^ (a_index % 2))
            in_b = not ((scan < b_endpoints[b_index]) ^ (b_index % 2))
            in_res = op(in_a, in_b)

            if in_res ^ (len(res) % 2):
                res += [scan]
            if scan == a_endpoints[a_index]:
                a_index += 1
            if scan == b_endpoints[b_index]:
                b_index += 1

            scan = min(a_endpoints[a_index], b_endpoints[b_index])

        # Reshape the resulting list of intervals to a Nx2 numpy array
        return np.asarray(res, dtype=np.int64).reshape((-1, 2))
