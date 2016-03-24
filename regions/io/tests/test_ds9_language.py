from __future__ import absolute_import, division, print_function, \
    unicode_literals
from numpy.testing import assert_allclose
from ..ds9_language import ds9_parser, region_list_to_objects, \
    objects_to_ds9_string
from astropy.utils.data import get_pkg_data_filename
from astropy.tests.helper import pytest


@pytest.mark.parametrize('filename',
                         ['data/ds9.fk5.reg',
                          'data/ds9.fk5.hms.reg',
                          'data/ds9.fk5.hms.strip.reg',
                          'data/ds9.fk5.strip.reg'])
def test_fk5(filename):
    filename = get_pkg_data_filename(filename)
    temp = ds9_parser(filename)
    regs = region_list_to_objects(temp)

    actual = objects_to_ds9_string(regs, coordsys='fk5', fmt='.2f', radunit='arcsec')

    # Use this to produce reference file for now
    #print(actual)
    #1/0

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
    temp = ds9_parser(filename)
    regs = region_list_to_objects(temp)

    actual = objects_to_ds9_string(regs, coordsys='galactic', fmt='.2f', radunit='arcsec')

    # Use this to produce reference file for now
    #print(actual)
    #1 / 0

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
    temp = ds9_parser(filename)
    regs = region_list_to_objects(temp)

    actual = objects_to_ds9_string(regs, coordsys='physical', fmt='.2f')

    # Use this to produce reference file for now
    # print(actual)
    # 1 / 0

    reference_file = get_pkg_data_filename('data/physical_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired

