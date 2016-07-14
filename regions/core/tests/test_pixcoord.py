# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_equal
from astropy.tests.helper import pytest
from ..pixcoord import PixCoord

try:
    import shapely

    HAS_SHAPELY = True
except:
    HAS_SHAPELY = False


def test_pixcoord_scalar():
    p = PixCoord(x=1, y=2)

    assert p.x == 1
    assert p.y == 2

    assert p.isscalar is True

    assert str(p) == 'PixCoord\nx : 1\ny : 2'
    assert repr(p) == 'PixCoord\nx : 1\ny : 2'

    with pytest.raises(TypeError):
        len(p)

    with pytest.raises(IndexError):
        p[0]


def test_pixcoord_array():
    p = PixCoord(x=[1, 2, 3], y=[11, 22, 33])

    assert_equal(p.x, [1, 2, 3])
    assert_equal(p.y, [11, 22, 33])

    assert p.isscalar is False

    assert str(p) == 'PixCoord\nx : [1 2 3]\ny : [11 22 33]'
    assert repr(p) == 'PixCoord\nx : [1 2 3]\ny : [11 22 33]'

    assert len(p) == 3

    # Test `__iter__` via assertions on the last element
    p2 = [_ for _ in p][-1]
    assert p2.x == 3
    assert p2.y == 33

    # Test `__getitem__
    p3 = p[-1]
    assert p3.isscalar
    assert p3.x == 3
    assert p3.y == 33

    p4 = p[1:]
    assert len(p4) == 2
    assert_equal(p4.x, [2, 3])
    assert_equal(p4.y, [22, 33])


@pytest.mark.skipif('not HAS_SHAPELY')
def test_pixcoord_shapely():
    from shapely.geometry.point import Point
    p = PixCoord(1, 2)
    s = p.to_shapely()
    assert isinstance(s, Point)
    assert s.x == 1
    assert s.y == 2
