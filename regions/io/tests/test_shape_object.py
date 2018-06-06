# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import pytest

from astropy.tests.helper import assert_quantity_allclose

from .. import DS9Parser, CRTFParser, to_shape_list


def test_shape_ds9():

    reg_str = "galactic\ncircle(42,43,3) # color=green"

    parser = DS9Parser(reg_str)
    shape1 = parser.shapes[0]
    region = parser.shapes.to_regions()
    shape2 = to_shape_list(region,'DS9','galactic')[0]

    assert shape1.format_type == shape2.format_type
    assert shape1.coordsys == shape2.coordsys
    assert shape1.region_type == shape2.region_type
    assert shape1.include == shape2.include
    assert shape1.coord == shape2.coord
    assert dict(shape1.meta) == dict(shape2.meta)
    assert shape1.composite == shape2.composite


def test_shape_crtf():

    reg_str = "circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, coord=J2000, symsize=2"

    parser = CRTFParser(reg_str)
    shape1 = parser.shapes[0]
    region = parser.shapes.to_regions()
    shape2 = to_shape_list(region, 'CRTF', 'fk5')[0]

    assert shape1.format_type == shape2.format_type
    assert shape1.coordsys == shape2.coordsys
    assert shape1.region_type == shape2.region_type
    assert shape1.include == shape2.include
    for x, y in zip(shape1.coord, shape2.coord):
        assert_quantity_allclose(x.to('deg'), y.to('deg'))
    # assert dict(shape1.meta) == dict(shape2.meta)
    assert shape1.composite == shape2.composite


def test_valid_shape():

    reg_str = "circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, coord=J2000, symsize=2"

    shape = CRTFParser(reg_str).shapes[0]
    with pytest.raises(ValueError) as err:
        shape.format_type = 'CAIO'

    assert "CAIO is not available as io" in str(err)

    shape = CRTFParser(reg_str).shapes[0]
    with pytest.raises(ValueError) as err:
        shape.region_type = 'box'

    assert "box is not a valid region type in this package" in str(err)

    shape = CRTFParser(reg_str).shapes[0]
    with pytest.raises(ValueError) as err:
        shape.coordsys = 'hello'

    assert "hello is not a valid coordinate reference frame in astropy" in str(err)
