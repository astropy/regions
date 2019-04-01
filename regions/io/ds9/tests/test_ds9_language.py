# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import distutils.version as vers

from numpy.testing import assert_allclose
import pytest

from astropy.utils.data import get_pkg_data_filename, get_pkg_data_filenames
import astropy.version as astrov
from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import catch_warnings, assert_quantity_allclose
from astropy.utils.exceptions import AstropyUserWarning
from astropy import units as u

from ....shapes.circle import CircleSkyRegion
from ..read import read_ds9, DS9Parser
from ..write import write_ds9, ds9_objects_to_string

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)


@pytest.mark.xfail(_ASTROPY_VERSION < _ASTROPY_MINVERSION,
                   reason='Some coordinates systems not available in older version of astropy')
def test_read():
    # Check that all test files including reference files are readable
    files = get_pkg_data_filenames('data')
    for f in files:
        with open(f) as f:
            DS9Parser(f.read(), errors='warn')


implemented_region_types = ('ellipse', 'circle', 'rectangle', 'polygon', 'point')


@pytest.mark.parametrize('filename',
                         ['data/ds9.fk5.reg',
                          'data/ds9.fk5.hms.reg',
                          'data/ds9.fk5.hms.strip.reg',
                          'data/ds9.fk5.strip.reg',
                          'data/ds9.galactic.reg',
                          'data/ds9.galactic.hms.reg',
                          'data/ds9.galactic.hms.strip.reg',
                          'data/ds9.galactic.strip.reg',
                          # TODO : data/ds9.physical.windows.reg contains different values -> Why?
                          'data/ds9.physical.reg',
                          'data/ds9.physical.strip.reg',
                         ])
def test_file(filename):
    filename = get_pkg_data_filename(filename)
    regs = read_ds9(filename, errors='warn')

    coordsys = os.path.basename(filename).split(".")[1]

    actual = ds9_objects_to_string(regs, coordsys=str(coordsys), fmt='.2f',
                                   radunit=None if coordsys=='physical' else 'arcsec').strip()

    reffile = get_pkg_data_filename('data/{coordsys}{strip}_reference.reg'
                                    .format(coordsys=coordsys,
                                            strip=("_strip"
                                                   if "strip" in filename
                                                   else "")))

    with open(reffile, 'r') as fh:
        desired = fh.read().strip()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split(" ")) for line in desired.split("\n")]
    actual_lines = [set(line.split(" ")) for line in actual.split("\n")]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


def test_ds9_objects_to_str():
    """
    Simple test case for ds9_objects_to_str.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    actual = ds9_objects_to_string([region], fmt='.4f')
    assert actual == expected


def test_ds9_string_to_objects():
    """
    Simple test case for ds9_string_to_objects.
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    parser = DS9Parser(ds9_str)
    regions = parser.shapes.to_regions()
    reg = regions[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)


def test_ds9_string_to_objects_whitespace():
    """
    Simple test case for ds9_string_to_objects.
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\n -circle(42.0000,43.0000,3.0000)\n'
    parser = DS9Parser(ds9_str)
    assert parser.shapes[0].include == False


def test_ds9_io(tmpdir):
    """
    Simple test case for write_ds9 and read_ds9.
    """
    center = SkyCoord(42, 43, unit='deg', frame='fk5')
    radius = Angle(3, 'deg')
    reg = CircleSkyRegion(center, radius)
    reg.meta['name'] = 'MyName'

    filename = os.path.join(str(tmpdir), 'ds9.reg')
    write_ds9([reg], filename, coordsys='fk5')
    reg = read_ds9(filename)[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)
    assert 'name' in reg.meta
    assert reg.meta['name'] == 'MyName'


def test_missing_region_warns():
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\nnotaregiontype(blah)'

    # this will warn on both the commented first line and the not_a_region line
    with catch_warnings(AstropyUserWarning) as ASWarn:
        parser = DS9Parser(ds9_str, errors='warn')

    assert len(parser.shapes) == 1
    assert len(ASWarn) == 1
    assert "Region type 'notaregiontype' was identified, but it is not one of the known region types." in str(ASWarn[0].message)


def test_global_parser():
    """
    Check that the global_parser does what's expected
    """

    # have to force "str" here because unicode_literals makes this string a
    # unicode string, even though in python2.7 we don't want it to be.  But in
    # py3, we don't want it to be bytes.
    global_test_str = str('global color=green dashlist=8 3 width=1'
                          ' font="helvetica 10 normal roman" select=1'
                          ' highlite=1 dash=0 fixed=0 edit=1 move=1'
                          ' delete=1 include=1 source=1')
    global_parser = DS9Parser(global_test_str)
    assert dict(global_parser.global_meta) == {'dash': '0', 'source': '1', 'move': '1',
                                  'font': 'helvetica 10 normal roman',
                                  'dashlist': '8 3', 'include': True,
                                  'highlite': '1', 'color': 'green',
                                  'select': '1',
                                  'fixed': '0', 'width': '1', 'edit': '1',
                                  'delete': '1'}


def test_ds9_color():
    """
    Color parsing test
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000) # color=green' \
              '\ncircle(43.0000,43.0000,3.0000) # color=orange\n'

    parser = DS9Parser(ds9_str)
    regions = parser.shapes

    assert regions[0].meta['color'] == 'green'
    assert regions[1].meta['color'] == 'orange'


def test_ds9_color_override_global():
    """
    Color parsing test in the presence of a global
    """
    global_test_str = str('global color=blue dashlist=8 3 width=1'
                          ' font="helvetica 10 normal roman" select=1'
                          ' highlite=1 dash=0 fixed=0 edit=1 move=1'
                          ' delete=1 include=1 source=1')

    ds9_str = '# Region file format: DS9 astropy/regions\n{global_str}\nfk5\n'
    reg1str = "circle(42.0000,43.0000,3.0000) # color=green"
    reg2str = "circle(43.0000,43.0000,3.0000) # color=orange"
    reg3str = "circle(43.0000,43.0000,3.0000)"
    ds9_str = ds9_str.format(global_str=global_test_str) + "\n".join([reg1str, reg2str, reg3str])
    parser = DS9Parser(ds9_str)
    regions = parser.shapes
    assert regions[0].meta['color'] == 'green'
    assert regions[1].meta['color'] == 'orange'
    assert regions[2].meta['color'] == 'blue'


def test_issue134_regression():
    regstr = 'galactic; circle(+0:14:26.064,+0:00:45.206,30.400")'
    parser = DS9Parser(regstr)
    regions = parser.shapes
    assert regions[0].to_region().radius.value == 30.4


def test_issue65_regression():
    regstr = 'J2000; circle 188.5557102 12.0314056 1" # color=red'
    parser = DS9Parser(regstr)
    regions = parser.shapes
    reg = regions[0].to_region()
    assert reg.center.ra.value == 188.5557102
    assert reg.center.dec.value == 12.0314056
    assert reg.radius.value == 1.0


def test_pixel_angle():
    """
    Checks whether angle in PixelRegions is a u.Quantity object.
    """
    reg_str = 'image\nbox(1.5,2,2,1,0)'

    reg = DS9Parser(reg_str).shapes.to_regions()[0]

    assert isinstance(reg.angle, u.quantity.Quantity)


def test_expicit_formatting_directives():
    """
    Checks whether every explicit formatting directive is supported
    """

    reg_str = 'image\ncircle(1.5, 2, 2)'
    valid_coord = DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'image\ncircle(1.5p, 2p, 2p)'
    coord = DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str_explicit = 'image\ncircle(1.5i, 2i, 2i)'
    coord = DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str = 'fk5\ncircle(1.5, 2, 2)'
    valid_coord = DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'fk5\ncircle(1.5d, 2d, 2d)'
    coord = DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str_explicit = 'fk5\ncircle(1.5r, 2r, 2r)'
    coord = DS9Parser(reg_str_explicit).shapes[0].coord
    for val in coord:
        assert val.unit == u.Unit('rad')

    reg_str = 'fk5\ncircle(1:20:30, 2:3:7, 2)'
    valid_coord = DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'fk5\ncircle(1h20m30s, 2d3m7s, 2d)'
    coord = DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(u.Quantity(coord), u.Quantity(valid_coord))

def test_text_metadata():
    """
    Regression test for issue #233: make sure that text metadata is parsed and
    stored appropriately.
    """

    reg_str = 'image\ncircle(1.5, 2, 2) # text={this_is_text}'
    parsed_reg = DS9Parser(reg_str)

    assert len(parsed_reg.shapes) == 1
    assert parsed_reg.shapes[0].meta['text'] == 'this_is_text'

    regs = parsed_reg.shapes.to_regions()

    assert regs[0].meta['label'] == 'this_is_text'
    assert regs[0].meta['text'] == 'this_is_text'
