# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the ds9 subpackage.
"""

import os

from numpy.testing import assert_allclose
import pytest

from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import catch_warnings, assert_quantity_allclose
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename, get_pkg_data_filenames
from astropy.utils.exceptions import AstropyUserWarning

from ....core import Regions
from ....shapes.circle import CircleSkyRegion
from ..core import DS9RegionParserWarning
from ..read import _DS9Parser


def test_read():
    # Check that all test files including reference files are readable
    files = get_pkg_data_filenames('data')
    # DS9RegionParserWarning from non-supported panda, [b/e]panda regions
    with catch_warnings(DS9RegionParserWarning):
        for file in files:
            with open(file) as fh:
                _DS9Parser(fh.read(), errors='warn')


# TODO: ds9.physical.windows.reg contains different values -> Why?
@pytest.mark.parametrize('filename',
                         ['data/ds9.fk5.reg',
                          'data/ds9.fk5.hms.reg',
                          'data/ds9.fk5.hms.strip.reg',
                          'data/ds9.fk5.strip.reg',
                          'data/ds9.galactic.reg',
                          'data/ds9.galactic.hms.reg',
                          'data/ds9.galactic.hms.strip.reg',
                          'data/ds9.galactic.strip.reg',
                          'data/ds9.physical.reg',
                          'data/ds9.physical.strip.reg'])
def test_file(filename):
    filename = get_pkg_data_filename(filename)

    # DS9RegionParserWarnings from skipped multi-annulus regions
    with catch_warnings(DS9RegionParserWarning):
        regs = Regions.read(filename, errors='warn', format='ds9')

    coordsys = os.path.basename(filename).split(".")[1]

    radunit = None if coordsys == 'physical' else 'arcsec'
    actual = regs.serialize(format='ds9', coordsys=str(coordsys), fmt='.2f',
                            radunit=radunit).strip()

    strip = '_strip' if 'strip' in filename else ''
    filepath = os.path.join('data', f'{coordsys}{strip}_reference.reg')
    reffile = get_pkg_data_filename(filepath)
    with open(reffile, 'r') as fh:
        desired = fh.read().strip()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split(' ')) for line in desired.split('\n')]
    actual_lines = [set(line.split(' ')) for line in actual.split('\n')]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


def test_ds9_serialize():
    """
    Simple test case for serialize.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = ('# Region file format: DS9 astropy/regions\nfk5\n'
                'circle(42.0000,43.0000,3.0000)\n')
    actual = region.serialize(format='ds9', fmt='.4f')
    assert actual == expected


def test_ds9_string_to_objects():
    """
    Simple test case for ds9_string_to_objects.
    """
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               'circle(42.0000,43.0000,3.0000)\n')
    regions = Regions.parse(ds9_str, format='ds9')
    reg = regions[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)


def test_ds9_string_to_objects_whitespace():
    """
    Simple test case for ds9_string_to_objects.
    """
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               ' -circle(42.0000,43.0000,3.0000)\n')
    region = Regions.parse(ds9_str, format='ds9')[0]
    assert not region.meta['include']


def test_ds9_io(tmpdir):
    """
    Simple test case for write and read.
    """
    center = SkyCoord(42, 43, unit='deg', frame='fk5')
    radius = Angle(3, 'deg')
    reg = CircleSkyRegion(center, radius)
    reg.meta['name'] = 'MyName'

    filename = os.path.join(str(tmpdir), 'ds9.reg')
    reg.write(filename, coordsys='fk5', format='ds9')
    reg2 = Regions.read(filename, format='ds9')[0]

    assert_allclose(reg2.center.ra.deg, 42)
    assert_allclose(reg2.center.dec.deg, 43)
    assert_allclose(reg2.radius.value, 3)
    assert 'name' in reg2.meta
    assert reg2.meta['name'] == 'MyName'


def test_missing_region_warns():
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               'circle(42.0000,43.0000,3.0000)\nnotaregiontype(blah)')

    # this will warn on both the commented first line and the
    # not_a_region line
    with catch_warnings(AstropyUserWarning) as ASWarn:
        regions = Regions.parse(ds9_str, format='ds9', errors='warn')

    assert len(regions) == 1
    assert len(ASWarn) == 1
    estr = ('Region type "notaregiontype" was found, but it is not one '
            'of the supported region types.')
    assert estr in str(ASWarn[0].message)


def test_global_parser():
    """
    Check that the global_parser does what's expected
    """
    global_test_str = ('global color=green dashlist=8 3 width=1 '
                       'font="helvetica 10 normal roman" select=1 '
                       'highlite=1 dash=0 fixed=0 edit=1 move=1 '
                       'delete=1 include=1 source=1')
    global_parser = _DS9Parser(global_test_str)

    exp = {'dash': '0', 'source': '1', 'move': '1',
           'font': 'helvetica 10 normal roman', 'dashlist': '8 3',
           'include': True, 'highlite': '1', 'color': 'green', 'select': '1',
           'fixed': '0', 'width': '1', 'edit': '1', 'delete': '1'}
    assert dict(global_parser.global_meta) == exp


def test_ds9_color():
    """
    Color parsing test
    """
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               'circle(42.0000,43.0000,3.0000) # color=green\n'
               'circle(43.0000,43.0000,3.0000) # color=orange\n')

    regions = Regions.parse(ds9_str, format='ds9')
    assert regions[0].visual['color'] == 'green'
    assert regions[1].visual['color'] == 'orange'


def test_ds9_color_override_global():
    """
    Color parsing test in the presence of a global
    """
    global_test_str = ('global color=blue dashlist=8 3 width=1 '
                       'font="helvetica 10 normal roman" select=1 '
                       'highlite=1 dash=0 fixed=0 edit=1 move=1 '
                       'delete=1 include=1 source=1\n')

    reg1str = 'circle(42.0000,43.0000,3.0000) # color=green'
    reg2str = 'circle(43.0000,43.0000,3.0000) # color=orange'
    reg3str = 'circle(43.0000,43.0000,3.0000)'

    global_str = (global_test_str + 'fk5\n'
                  + '\n'.join([reg1str, reg2str, reg3str]))
    ds9_str = f'# Region file format: DS9 astropy/regions\n{global_str}\n'
    regions = Regions.parse(ds9_str, format='ds9')
    assert regions[0].visual['color'] == 'green'
    assert regions[1].visual['color'] == 'orange'
    assert regions[2].visual['color'] == 'blue'


def test_issue134_regression():
    regstr = 'galactic; circle(+0:14:26.064,+0:00:45.206,30.400")'
    regions = Regions.parse(regstr, format='ds9')
    assert regions[0].radius.value == 30.4


def test_issue65_regression():
    regstr = 'J2000; circle 188.5557102 12.0314056 1" # color=red'
    region = Regions.parse(regstr, format='ds9')[0]
    assert region.center.ra.value == 188.5557102
    assert region.center.dec.value == 12.0314056
    assert region.radius.value == 1.0


def test_frame_uppercase():
    """
    Regression test for issue #236 (PR #237)
    """
    regstr = ('GALACTIC\ncircle(188.5557102,12.0314056,7.1245) # '
              'text={This message has both ' 'a " and : in it} textangle=30')
    region = Regions.parse(regstr, format='ds9')[0]
    assert region.center.l.value == 188.5557102
    assert region.center.b.value == 12.0314056
    assert region.radius.value == 7.1245


def test_pixel_angle():
    """
    Check whether angle in PixelRegions is a u.Quantity object.
    """
    reg_str = 'image\nbox(1.5,2,2,1,0)'
    region = Regions.parse(reg_str, format='ds9')[0]
    assert isinstance(region.angle, u.quantity.Quantity)


def test_expicit_formatting_directives():
    """
    Check whether every explicit formatting directive is supported.
    """
    reg_str = 'image\ncircle(1.5, 2, 2)'
    valid_coord = _DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'image\ncircle(1.5p, 2p, 2p)'
    coord = _DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str_explicit = 'image\ncircle(1.5i, 2i, 2i)'
    coord = _DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str = 'fk5\ncircle(1.5, 2, 2)'
    valid_coord = _DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'fk5\ncircle(1.5d, 2d, 2d)'
    coord = _DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(coord, valid_coord)

    reg_str_explicit = 'fk5\ncircle(1.5r, 2r, 2r)'
    coord = _DS9Parser(reg_str_explicit).shapes[0].coord
    for val in coord:
        assert val.unit == u.Unit('rad')

    reg_str = 'fk5\ncircle(1:20:30, 2:3:7, 2)'
    valid_coord = _DS9Parser(reg_str).shapes[0].coord

    reg_str_explicit = 'fk5\ncircle(1h20m30s, 2d3m7s, 2d)'
    coord = _DS9Parser(reg_str_explicit).shapes[0].coord
    assert_quantity_allclose(u.Quantity(coord), u.Quantity(valid_coord))


def test_text_metadata():
    """
    Regression test for issue #233: make sure that text metadata is
    parsed and stored appropriately.
    """
    reg_str = 'image\ncircle(1.5, 2, 2) # text={this_is_text}'
    regions = Regions.parse(reg_str, format='ds9')

    assert len(regions) == 1
    assert regions[0].meta['text'] == 'this_is_text'
    assert regions[0].meta['label'] == 'this_is_text'
    assert regions[0].meta['text'] == 'this_is_text'


def test_angle_serialization():
    """
    Regression test for issue #223 to ensure Angle arcsec inputs are
    correctly converted to degrees.
    """
    reg = Regions([CircleSkyRegion(SkyCoord(10,20, unit='deg'),
                                   Angle(1, 'arcsec'))])
    regstr = reg.serialize(format='ds9')
    expected = ('# Region file format: DS9 astropy/regions\nfk5\n'
                'circle(10.000009,20.000002,0.000278)\n')
    assert regstr == expected


def test_semicolon():
    """
    Regression test for issue #238 to allow semicolons in the text
    field.
    """
    regstr = ('galactic\n'
              'circle(0.003, 0.1, 206.696") # text={S17; test} color=red\n'
              'circle(0.1, -0.5, 360.148") # text={S19} color=red')

    reg = Regions.parse(regstr, format='ds9')
    assert len(reg) == 2
