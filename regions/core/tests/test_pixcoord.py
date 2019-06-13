# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from numpy.testing import assert_equal, assert_allclose
import pytest
from astropy import units as u
from ..._utils.examples import make_example_dataset
from ..pixcoord import PixCoord


@pytest.fixture(scope='session')
def wcs():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs


def test_pixcoord_copy_scalar():
    p = PixCoord(x=1, y=2)
    p2 = p.copy()
    p.x = 99
    p.y = 99
    assert p2.xy == (1, 2)


def test_pixcoord_copy_array():
    p = PixCoord(x=[1, 2], y=[3, 4])
    p2 = p.copy()
    p.x[0] = 99
    p.y[0] = 99
    assert_equal(p2.x, [1, 2])
    assert_equal(p2.y, [3, 4])


def test_pixcoord_basic_dimension():
    with pytest.raises(ValueError):
        PixCoord(np.array([1, 2]), [3, 4, 5, 6])


def test_pixcoord_basics_scalar():
    p = PixCoord(x=1, y=2)
    p1 = PixCoord(x=np.array(1), y=2)
    p2 = PixCoord(x=np.array(1), y=np.array(2))

    assert p == p1
    assert p1 == p2

    assert p.x == 1
    assert p.y == 2

    assert p.isscalar
    assert p1.isscalar
    assert p2.isscalar

    assert str(p) == 'PixCoord(x=1, y=2)'
    assert repr(p) == 'PixCoord(x=1, y=2)'

    with pytest.raises(TypeError):
        len(p)

    with pytest.raises(IndexError):
        p[0]


def test_pixcoord_basics_array_1d():
    p = PixCoord(x=[1, 2, 3], y=[11, 22, 33])
    p1 = PixCoord(x=1, y=[1, 2, 3])

    assert_equal(p.x, [1, 2, 3])
    assert_equal(p.y, [11, 22, 33])

    assert_equal(p1.x, [1, 1, 1])
    assert_equal(p1.y, [1, 2, 3])

    assert not p.isscalar
    assert not p1.isscalar

    assert str(p) == 'PixCoord(x=[1 2 3], y=[11 22 33])'
    assert str(p1) == 'PixCoord(x=[1 1 1], y=[1 2 3])'

    assert repr(p) == 'PixCoord(x=[1 2 3], y=[11 22 33])'
    assert repr(p1) == 'PixCoord(x=[1 1 1], y=[1 2 3])'

    assert len(p) == 3
    assert len(p1) == 3

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


def test_pixcoord_basics_array_2d():
    p = PixCoord(
        [[1, 2, 3], [4, 5, 6]],
        [[11, 12, 13], [14, 15, 16]],
    )

    assert_equal(p.x, [[1, 2, 3], [4, 5, 6]])
    assert_equal(p.y, [[11, 12, 13], [14, 15, 16]])

    assert p.isscalar is False

    # TODO: this repr with the newline isn't nice ... improve it somehow?
    expected = ('PixCoord(x=[[1 2 3]\n [4 5 6]], '
                'y=[[11 12 13]\n [14 15 16]])')
    assert str(p) == expected
    assert repr(p) == expected

    assert len(p) == 2

    # TODO: does `__iter__` and `__getitem__` make sense for 2d array case?
    # (see tests above)


def test_pixcoord_to_sky_scalar(wcs):
    p = PixCoord(x=18, y=9)

    s = p.to_sky(wcs=wcs)
    assert s.name == 'galactic'
    assert_allclose(s.data.lon.deg, 349.88093995435634)
    assert_allclose(s.data.lat.deg, 10.003028030623508)

    p2 = PixCoord.from_sky(skycoord=s, wcs=wcs)
    assert isinstance(p2.x, float)
    assert p2.isscalar

    assert p == p2


def test_pixcoord_to_sky_array_1d(wcs):
    p = PixCoord(x=[17, 18], y=[8, 9])

    s = p.to_sky(wcs=wcs)
    assert s.name == 'galactic'
    assert_allclose(s.data.lon.deg, [0, 349.88094])
    assert_allclose(s.data.lat.deg, [0, 10.003028])

    p2 = PixCoord.from_sky(skycoord=s, wcs=wcs)
    assert isinstance(p2.x, np.ndarray)
    assert p2.x.shape == (2,)
    assert not p2.isscalar

    assert p == p2


def test_pixcoord_to_sky_array_2d(wcs):
    p1 = PixCoord(x=[[17, 17, 17], [18, 18, 18]], y=[[8, 8, 8], [9, 9, 9]])
    p = PixCoord(x=[[17, 18]], y=[[8, 9]])

    s = p.to_sky(wcs=wcs)
    assert s.name == 'galactic'

    s1 = p1.to_sky(wcs=wcs)
    assert s.name == 'galactic'
    assert_allclose(s1.data.lon.deg, [[0, 0, 0], [349.88094, 349.88094, 349.88094]])
    assert_allclose(s1.data.lat.deg, [[0, 0, 0], [10.003028, 10.003028, 10.003028]])
    assert_allclose(s.data.lon.deg, [[0, 349.88094]])
    assert_allclose(s.data.lat.deg, [[0, 10.003028]])

    p2 = PixCoord.from_sky(skycoord=s, wcs=wcs)
    assert isinstance(p2.x, np.ndarray)
    assert p2.x.shape == (1, 2)
    assert not p2.isscalar

    assert p == p2


def test_pixcoord_separation_scalar():
    p1 = PixCoord(x=1, y=2)
    p2 = PixCoord(x=4, y=6)
    sep = p1.separation(p2)
    assert_allclose(sep, 5)


def test_pixcoord_separation_array_1d():
    p1 = PixCoord(x=[1, 1], y=[2, 2])
    p2 = PixCoord(x=[4, 4], y=[6, 6])
    sep = p1.separation(p2)
    assert_allclose(sep, [5, 5])


def test_pixcoord_separation_array_2d():
    p1 = PixCoord(x=[[1, 1]], y=[[2, 2]])
    p2 = PixCoord(x=[[4, 4]], y=[[6, 6]])
    sep = p1.separation(p2)
    assert isinstance(sep, np.ndarray)
    assert sep.shape == (1, 2)
    assert_allclose(sep, [[5, 5]])


def test_equality():
    a = np.array([1, 2])
    b = PixCoord(a[0], a[1])
    c = PixCoord(a[0] + 0.0000001, a[1])

    assert not b == a
    assert b == b
    assert b == c

    a = PixCoord(
        [[1, 2, 3], [4, 5, 6]],
        [[11, 12, 13], [14, 15, 16]],
    )

    b = PixCoord(
        [[1, 2, 3], [4, 5, 6]],
        [[11.0000002, 12, 13], [14, 15, 16]],
    )

    assert a == b
    assert a == a


def test_pixcoord_xy():
    a = np.array([1, 2])
    pc = PixCoord(a[0], a[1])

    assert pc.xy[0] == pc.x
    assert pc.xy[1] == pc.y


def test_pixcoord_rotate_scalar():
    point = PixCoord(3, 4)
    center = PixCoord(2, 3)
    point = point.rotate(center, 90 * u.deg)
    assert_allclose(point.xy, (1, 4))


def test_pixcoord_rotate_array():
    point = PixCoord([3, 3], [4, 4])
    center = PixCoord(2, 3)
    point = point.rotate(center, 90 * u.deg)
    assert_allclose(point.xy, ([1, 1], [4, 4]))
