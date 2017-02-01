# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import distutils.version as vers

from numpy.testing import assert_allclose

from astropy.utils.data import get_pkg_data_filename, get_pkg_data_filenames
from astropy.tests.helper import pytest
import astropy.version as astrov
from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import catch_warnings
from astropy.utils.exceptions import AstropyUserWarning

from ...shapes.circle import CircleSkyRegion
from ..read_ds9 import read_ds9, ds9_string_to_objects, global_parser
from ..write_ds9 import write_ds9, ds9_objects_to_string

_ASTROPY_MINVERSION = vers.LooseVersion('1.1')
_ASTROPY_VERSION = vers.LooseVersion(astrov.version)


@pytest.mark.xfail(_ASTROPY_VERSION < _ASTROPY_MINVERSION,
                   reason='Some coordinates systems not available in older version of astropy')
def test_read():
    # Check that all test files including reference files are readable
    files = get_pkg_data_filenames('data')
    for f in files:
        read_ds9(f, errors='warn')


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

    reference_file = get_pkg_data_filename('data/fk5_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired


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

    reference_file = get_pkg_data_filename('data/galactic_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired


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

    reference_file = get_pkg_data_filename('data/physical_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired


def test_ds9_objects_to_str():
    """Simple test case for ds9_objects_to_str.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    actual = ds9_objects_to_string([region])
    assert actual == expected


def test_ds9_string_to_objects():
    """Simple test case for ds9_string_to_objects
    """
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    regions = ds9_string_to_objects(ds9_str)
    reg = regions[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)


def test_ds9_io(tmpdir):
    """Simple test case for write_ds9 and read_ds9.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    reg = CircleSkyRegion(center, radius)

    filename = os.path.join(str(tmpdir), 'ds9.reg')
    write_ds9([reg], filename)
    reg = read_ds9(filename)[0]

    assert_allclose(reg.center.ra.deg, 42)
    assert_allclose(reg.center.dec.deg, 43)
    assert_allclose(reg.radius.value, 3)


def test_missing_region_warns():
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\nnotaregiontype(blah)'

    # this will warn on both the commented first line and the not_a_region line
    with catch_warnings(AstropyUserWarning) as ASWarn:
        regions = ds9_string_to_objects(ds9_str, errors='warn')

    assert len(regions) == 1
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
    global_parsed = global_parser(global_test_str)
    assert global_parsed[0] == 'global'
    assert global_parsed[1] == {'dash': '0', 'source': '1', 'move': '1',
                                'font': '"helvetica 10 normal roman" ',
                                'dashlist': '8 3 ', 'include': '1',
                                'highlite': '1', 'color': 'green',
                                'select': '1',
                                'fixed': '0', 'width': '1', 'edit': '1',
                                'delete': '1'}

def test_ds9_color():
    """Color parsing test"""
    ds9_str = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000) # color=green\ncircle(43.0000,43.0000,3.0000) # color=orange\n'
    regions = ds9_string_to_objects(ds9_str)

    assert regions[0].visual['color'] == 'green'
    assert regions[1].visual['color'] == 'orange'

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
    regions = ds9_string_to_objects(ds9_str)

    assert regions[0].visual['color'] == 'green'
    assert regions[1].visual['color'] == 'orange'
    assert regions[2].visual['color'] == 'blue'
