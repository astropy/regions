# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy.units.quantity import Quantity

from .. import DS9Parser, CRTFParser, to_shape_list


def test_shape_ds9():

    reg_str = "galactic\ncircle(42,43,3) # color=green"

    parser = DS9Parser(reg_str)
    shape1 = parser.shapes[0]
    region = parser.shapes.to_regions()
    shape2 = to_shape_list(region, 'galactic')[0]

    assert shape1.coordsys == shape2.coordsys
    assert shape1.region_type == shape2.region_type
    assert shape1.include == shape2.include
    assert shape1.coord == shape2.coord
    assert dict(shape1.meta) == dict(shape2.meta)
    assert shape1.composite == shape2.composite

    # Checks for pixel offset since DS9 has a origin=1 pixel system
    reg_str = 'image\nbox(1, 2, 3, 4, 5)'

    parser = DS9Parser(reg_str)
    shape = parser.shapes[0]

    assert shape.coord[0] == 0
    assert shape.coord[1] == 1
    assert_quantity_allclose(shape.coord[2:-1], [Quantity("3"), Quantity("4")])
    assert_quantity_allclose(shape.coord[-1], Quantity("5deg"))


def test_shape_crtf():

    reg_str = "circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, coord=J2000, symsize=2"

    parser = CRTFParser(reg_str)
    shape1 = parser.shapes[0]
    region = parser.shapes.to_regions()
    shape2 = to_shape_list(region, 'fk5')[0]

    assert shape1.coordsys == shape2.coordsys
    assert shape1.region_type == shape2.region_type
    assert shape1.include == shape2.include
    assert_quantity_allclose(shape1.coord, shape2.coord)
    # assert dict(shape1.meta) == dict(shape2.meta)
    assert shape1.composite == shape2.composite

    # Pixel origin is assumed to be 0 origin. Needs to be checked though.
    reg_str = "circle[[1pix, 2pix], 3pix], linewidth=2, coord=image, symsize=2"

    parser = CRTFParser(reg_str)
    shape = parser.shapes[0]
    assert_quantity_allclose(shape.coord, [Quantity("1"), Quantity("2"), Quantity("3")])


def test_valid_shape():

    reg_str = "circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, coord=J2000, symsize=2"

    shape = CRTFParser(reg_str).shapes[0]

    with pytest.raises(ValueError) as err:
        shape.region_type = 'box'

    assert "'box' is not a valid region type in this package" in str(err)

    shape = CRTFParser(reg_str).shapes[0]
    with pytest.raises(ValueError) as err:
        shape.coordsys = 'hello'

    assert "'hello' is not a valid coordinate reference frame in astropy" in str(err)


def test_valid_ellipse_ds9():

    reg_str1 = "image\nellipse (1, 2, 3, 4, 5)"
    reg_str2 = "fk5\nellipse (1, 2, 3, 4, 5)"
    reg_str3 = "image\nellipse (1, 2, 3, 4, 5, 6, 7)"  # Elliptical Annulus

    shape1 = DS9Parser(reg_str1, 'warn').shapes[0]
    shape2 = DS9Parser(reg_str2, 'warn').shapes[0]
    shape3 = DS9Parser(reg_str3, 'warn').shapes[0]

    assert_quantity_allclose(shape1.coord[:2], [0, 1])
    assert_quantity_allclose(shape2.coord[:2], [Quantity("1deg"), Quantity("2deg")])
    assert_quantity_allclose(shape3.coord[:2], [0, 1])

    assert_quantity_allclose(shape1.coord[2:-1], [6, 8])
    assert_quantity_allclose(shape2.coord[2:-1], [Quantity("6deg"), Quantity("8deg")])
    assert_quantity_allclose(shape3.coord[2:-1], [6, 8, 10, 12])

    assert_quantity_allclose(shape1.coord[-1], Quantity("5deg"))
    assert_quantity_allclose(shape2.coord[-1], Quantity("5deg"))
    assert_quantity_allclose(shape3.coord[-1], Quantity("7deg"))


def test_valid_ellipse_crtf():

    reg_str1 = "ellipse[[1pix, 2pix], [3pix, 4pix], 5deg]"
    reg_str2 = "ellipse[[1deg, 2deg], [3deg, 4deg], 5deg], coord=fk5"

    shape1 = CRTFParser(reg_str1, 'warn').shapes[0]
    shape2 = CRTFParser(reg_str2, 'warn').shapes[0]

    assert_quantity_allclose(shape1.coord[:2], [1, 2])
    assert_quantity_allclose(shape2.coord[:2], [Quantity("1deg"), Quantity("2deg")])

    assert_quantity_allclose(shape1.coord[2:-1], [6, 8])
    assert_quantity_allclose(shape2.coord[2:-1], [Quantity("6deg"), Quantity("8deg")])

    assert_quantity_allclose(shape1.coord[-1], Quantity("5deg"))
    assert_quantity_allclose(shape2.coord[-1], Quantity("5deg"))
