# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_equal, assert_allclose
from astropy.tests.helper import pytest
from ...utils.examples import make_example_dataset
from ..pixcoord import PixCoord

try:
    import shapely

    HAS_SHAPELY = True
except:
    HAS_SHAPELY = False


@pytest.fixture(scope='session')
def wcs():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs


def test_pixcoord_scalar_basics():
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


def test_pixcoord_array_basics():
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


def test_pixcoord_scalar_sky(wcs):
    p = PixCoord(x=18, y=9)

    s = p.to_sky(wcs=wcs)
    assert s.name == 'galactic'
    assert_allclose(s.data.lon.deg, 349.88093995435634)
    assert_allclose(s.data.lat.deg, 10.003028030623508)

    p2 = PixCoord.from_sky(skycoord=s, wcs=wcs)
    assert_allclose(p2.x, 18)
    assert_allclose(p2.y, 9)


def test_pixcoord_array_sky(wcs):
    p = PixCoord(x=[17, 18], y=[8, 9])

    s = p.to_sky(wcs=wcs)
    assert s.name == 'galactic'
    assert_allclose(s.data.lon.deg, [0, 349.88094])
    assert_allclose(s.data.lat.deg, [0, 10.003028])

    p2 = PixCoord.from_sky(skycoord=s, wcs=wcs)
    assert_allclose(p2.x, [17, 18])
    assert_allclose(p2.y, [8, 9])


def test_pixcoord_separation():
    # check scalar
    p1 = PixCoord(x=1, y=2)
    p2 = PixCoord(x=4, y=6)
    sep = p1.separation(p2)
    assert_allclose(sep, 5)

    # check array
    p1 = PixCoord(x=[1, 1], y=[2, 2])
    p2 = PixCoord(x=[4, 4], y=[6, 6])
    sep = p1.separation(p2)
    assert_allclose(sep, [5, 5])


@pytest.mark.skipif('not HAS_SHAPELY')
def test_pixcoord_shapely():
    from shapely.geometry.point import Point
    p = PixCoord(x=1, y=2)
    s = p.to_shapely()
    assert isinstance(s, Point)
    assert s.x == 1
    assert s.y == 2

    p2 = PixCoord.from_shapely(point=s)
    assert p2.x == 1
    assert p2.y == 2

    p = PixCoord(x=[1, 2, 3], y=[11, 22, 33])
    with pytest.raises(TypeError):
        p.to_shapely()
