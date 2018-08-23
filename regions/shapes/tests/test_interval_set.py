import pytest
import numpy as np
from ..moc import IntervalSet


@pytest.fixture()
def isets():
    a = IntervalSet(np.asarray([(49, 73), (53, 54), (33, 63), (65, 80), (51, 80), (100, 126), (38, 68), (61, 72), (74, 102), (27, 43)]))
    b = IntervalSet(np.asarray([(17, 26), (17, 41), (12, 31), (32, 61), (68, 90), (77, 105), (18, 27), (12, 35), (9, 37), (87, 97)]))
    return dict(a=a, b=b)


def test_interval_set_consistency(isets):
    assert isets['a'] == IntervalSet(np.asarray([(27, 126)]))
    assert isets['b'] == IntervalSet(np.asarray([(9, 61), (68, 105)]))


def test_interval_set_union(isets):
    assert isets['a'].union(isets['b']) == IntervalSet(np.array([(9, 126)]))
    assert isets['a'].union(IntervalSet()) == IntervalSet(np.array([(27, 126)]))
    assert IntervalSet().union(isets['a']) == IntervalSet(np.array([(27, 126)]))


def test_interval_set_intersection(isets):
    assert isets['a'].intersection(isets['b']) == IntervalSet(np.asarray([(27, 61), (68, 105)]))
    assert isets['a'].intersection(IntervalSet()) == IntervalSet()
    assert IntervalSet().intersection(isets['a']) == IntervalSet()


def test_interval_set_difference(isets):
    assert isets['a'].difference(isets['b']) == IntervalSet(np.asarray([(61, 68), (105, 126)]))
    assert isets['b'].difference(isets['a']) == IntervalSet(np.asarray([(9, 27)]))
    assert IntervalSet().difference(isets['a']) == IntervalSet()
    assert isets['a'].difference(IntervalSet()) == isets['a']


@pytest.fixture()
def isets2():
    nested1 = IntervalSet(np.array([[0, 1]]))
    nuniq1 = IntervalSet(np.array([[4*4**29, 4*4**29 + 1]]))
    nested2 = IntervalSet(np.array([[7, 76]]))
    nuniq2 = IntervalSet(np.array([[1 + 4*4**27, 4 + 4*4**27],
                                   [2 + 4*4**28, 4 + 4*4**28],
                                   [16 + 4*4**28, 19 + 4*4**28],
                                   [7 + 4*4**29,  8 + 4*4**29]]))
    return dict(nested1=nested1, nuniq1=nuniq1,
                nested2=nested2, nuniq2=nuniq2)


def test_to_nuinq_interval_set(isets2):
    assert IntervalSet.to_nuniq_interval_set(isets2['nested1']) == isets2['nuniq1']
    assert IntervalSet.to_nuniq_interval_set(isets2['nested2']) == isets2['nuniq2']
    # empty nested interval set
    assert IntervalSet.to_nuniq_interval_set(IntervalSet()) == IntervalSet()


def test_from_nuinq_interval_set(isets2):
    assert IntervalSet.from_nuniq_interval_set(isets2['nuniq1']) == isets2['nested1']
    assert IntervalSet.from_nuniq_interval_set(isets2['nuniq2']) == isets2['nested2']
    # empty nuniq interval set
    assert IntervalSet.from_nuniq_interval_set(IntervalSet()) == IntervalSet()


def test_from_to_interval_set(isets2):
    assert IntervalSet.from_nuniq_interval_set(
        IntervalSet.to_nuniq_interval_set(isets2['nested1'])
    ) == isets2['nested1']


def test_interval_set_min(isets):
    assert isets['a'].min == 27
    assert isets['b'].min == 9
    assert isets['a'].union(isets['b']).min == 9


def test_interval_set_max(isets):
    assert isets['a'].max == 126
    assert isets['b'].max == 105
    assert isets['a'].union(isets['b']).max == 126


def test_repr_interval_set(isets):
    assert repr(isets['a']) == "[[ 27 126]]"
    assert repr(isets['b']) == "[[  9  61]\n" \
                               " [ 68 105]]"
