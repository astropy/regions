# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy.tests.helper import pytest
from ..pixcoord import PixCoord

try:
    import shapely

    HAS_SHAPELY = True
except:
    HAS_SHAPELY = False


def test_pixcoord_basics():
    p = PixCoord(1, 2)

    assert p.x == 1
    assert p.y == 2

    # Test __iter__
    assert list(p) == [1, 2]

    # assert p.isscalar is True


@pytest.mark.skipif('not HAS_SHAPELY')
def test_pixcoord_shapely():
    from shapely.geometry.point import Point
    p = PixCoord(1, 2)
    s = p.to_shapely()
    assert isinstance(s, Point)
    assert s.x == 1
    assert s.y == 2
