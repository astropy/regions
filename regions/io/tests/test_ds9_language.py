# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import distutils.version as v
from astropy.utils.data import get_pkg_data_filename, get_pkg_data_filenames
from astropy.tests.helper import pytest
import astropy.version as astrov
from astropy.coordinates import Angle, SkyCoord
from ...shapes.circle import CircleSkyRegion
from ..read_ds9 import read_ds9
from ..write_ds9 import objects_to_ds9_string

_ASTROPY_MINVERSION = v.LooseVersion('1.1')
_ASTROPY_VERSION = v.LooseVersion(astrov.version)


@pytest.mark.xfail(_ASTROPY_VERSION < _ASTROPY_MINVERSION,
                   reason='Some coordinates systems not available in older version of astropy')
def test_read():
    # Check that all test files including reference files are readable
    files = get_pkg_data_filenames('data')
    for f in files:
        read_ds9(f)


@pytest.mark.parametrize('filename',
                         ['data/ds9.fk5.reg',
                          'data/ds9.fk5.hms.reg',
                          'data/ds9.fk5.hms.strip.reg',
                          'data/ds9.fk5.strip.reg'])
def test_fk5(filename):
    filename = get_pkg_data_filename(filename)
    regs = read_ds9(filename)

    actual = objects_to_ds9_string(regs, coordsys='fk5', fmt='.2f', radunit='arcsec')

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
    regs = read_ds9(filename)

    actual = objects_to_ds9_string(regs, coordsys='galactic', fmt='.2f', radunit='arcsec')

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
    regs = read_ds9(filename)

    actual = objects_to_ds9_string(regs, coordsys='physical', fmt='.2f')

    # Use this to produce reference file for now
    # print(actual)
    # 1 / 0

    reference_file = get_pkg_data_filename('data/physical_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired


def test_ds9_circle():
    """Test that circle to ds9 string works.

    Regression test for https://github.com/astropy/regions/issues/41
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = '# Region file format: DS9 astropy/regions\nfk5\ncircle(42.0000,43.0000,3.0000)\n'
    actual = objects_to_ds9_string([region])
    assert actual == expected
