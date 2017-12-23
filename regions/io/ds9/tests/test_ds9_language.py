# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import distutils.version as vers

from numpy.testing import assert_allclose

from astropy.utils.data import get_pkg_data_filename, get_pkg_data_filenames
import pytest
import astropy.version as astrov
from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import catch_warnings
from astropy.utils.exceptions import AstropyUserWarning

from ....shapes.circle import CircleSkyRegion
from .. import read_ds9, write_ds9, ds9_objects_to_string, DS9Parser

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)


@pytest.mark.xfail(_ASTROPY_VERSION < _ASTROPY_MINVERSION,
                   reason='Some coordinates systems not available in older version of astropy')
def test_read():
    # Check that all test files including reference files are readable
    files = get_pkg_data_filenames('data')
    for f in files:
        read_ds9(f, errors='warn')


implemented_region_types = ('ellipse', 'circle', 'rectangle', 'polygon', 'point')

@pytest.mark.parametrize('filename',
                         ['data/ds9.fk5.reg',
                          'data/ds9.fk5.hms.reg',
                          'data/ds9.fk5.hms.strip.reg',
                          'data/ds9.fk5.strip.reg'])
def test_fk5(filename):
    filename = get_pkg_data_filename(filename)
    regs = read_ds9(filename, errors='warn')

    actual = ds9_objects_to_string(regs, coordsys='fk5', fmt='.2f', radunit='arcsec')

    # Use this to produce reference file for now
    # print(actual)
    # 1/0

    if 'strip' in filename:
        reference_file = get_pkg_data_filename('data/fk5_strip_reference.reg')
    else:
        # preserves metadata
        reference_file = get_pkg_data_filename('data/fk5_reference.reg')

    with open(reference_file, 'r') as fh:
        desired = fh.read()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split()) for line in desired.split("\n")]
    actual_lines = [set(line.split()) for line in actual.split("\n")]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


@pytest.mark.parametrize('filename',
                         ['data/ds9.galactic.reg',
                          'data/ds9.galactic.hms.reg',
                          'data/ds9.galactic.hms.strip.reg',
                          'data/ds9.galactic.strip.reg'])
def test_galactic(filename):
    filename = get_pkg_data_filename(filename)
    regs = read_ds9(filename, errors='warn')

    actual = ds9_objects_to_string(regs, coordsys='galactic', fmt='.2f', radunit='arcsec')

    # Use this to produce reference file for now
    # print(actual)
    # 1 / 0

    if 'strip' in filename:
        reference_file = get_pkg_data_filename('data/galactic_strip_reference.reg')
    else:
        reference_file = get_pkg_data_filename('data/galactic_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split()) for line in desired.split("\n")]
    actual_lines = [set(line.split()) for line in actual.split("\n")]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


# Todo : data/ds9.physical.windows.reg contains different values -> Why?
@pytest.mark.parametrize('filename',
                         ['data/ds9.physical.reg',
                          'data/ds9.physical.strip.reg'])
def test_physical(filename):
    filename = get_pkg_data_filename(filename)
    regs = read_ds9(filename, errors='warn')

    actual = ds9_objects_to_string(regs, coordsys='physical', fmt='.2f')

    # Use this to produce reference file for now
    # print(actual)
    # 1 / 0

    if 'strip' in filename:
        reference_file = get_pkg_data_filename('data/physical_strip_reference.reg')
    else:
        reference_file = get_pkg_data_filename('data/physical_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split()) for line in desired.split("\n")]
    actual_lines = [set(line.split()) for line in actual.split("\n")]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


def test_ds9_objects_to_str():
    """Simple test case for ds9_objects_to_str.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    actual = ds9_objects_to_string([region], fmt='.4f')
    assert actual == expected


def test_ds9_string_to_objects():
    """Simple test case for ds9_string_to_objects
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    parser = DS9Parser(ds9_str)
    parser.run()
    regions = parser.shapes.to_region()
    reg = regions[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)

def test_ds9_string_to_objects_whitespace():
    """Simple test case for ds9_string_to_objects
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\n -circle(42.0000,43.0000,3.0000)\n'
    parser = DS9Parser(ds9_str)
    parser.run()
    assert parser.shapes[0].include == False

def test_ds9_io(tmpdir):
    """Simple test case for write_ds9 and read_ds9.
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
    assert "Region type 'notaregiontype'" in str(ASWarn[0].message)


def test_global_parser():
    """ Check that the global_parser does what's expected """
    # have to force "str" here because unicode_literals makes this string a
    # unicode string, even though in python2.7 we don't want it to be.  But in
    # py3, we don't want it to be bytes.
    global_test_str = str('global color=green dashlist=8 3 width=1'
                          ' font="helvetica 10 normal roman" select=1'
                          ' highlite=1 dash=0 fixed=0 edit=1 move=1'
                          ' delete=1 include=1 source=1')
    global_parser = DS9Parser(global_test_str)
    global_parser.run()
    global_parser.global_meta == {'dash': '0', 'source': '1', 'move': '1',
                                  'font': '"helvetica 10 normal roman" ',
                                  'dashlist': '8 3 ', 'include': '1',
                                  'highlite': '1', 'color': 'green',
                                  'select': '1',
                                  'fixed': '0', 'width': '1', 'edit': '1',
                                  'delete': '1'}

def test_ds9_color():
    """Color parsing test"""
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000) # color=green\ncircle(43.0000,43.0000,3.0000) # color=orange\n'

    parser = DS9Parser(ds9_str)
    parser.run()
    regions = parser.shapes

    assert regions[0].meta['color'] == 'green'
    assert regions[1].meta['color'] == 'orange'

def test_ds9_color_override_global():
    """Color parsing test in the presence of a global"""
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
    parser.run()
    regions = parser.shapes
    assert regions[0].meta['color'] == 'green'
    assert regions[1].meta['color'] == 'orange'
    assert regions[2].meta['color'] == 'blue'

def test_issue134_regression():
    regstr = 'galactic; circle(+0:14:26.064,+0:00:45.206,30.400")'
    parser = DS9Parser(regstr)
    parser.run()
    regions = parser.shapes
    assert regions[0].to_region().radius.value == 30.4
