from __future__ import absolute_import, division, print_function, unicode_literals

from numpy.testing import assert_allclose

from ..ds9_language import ds9_parser, region_list_to_objects, objects_to_ds9_string
from astropy.utils.data import get_pkg_data_filename


def test_fk5(tmpdir):
    filename = get_pkg_data_filename('data/ds9.fk5.reg')
    temp = ds9_parser(filename)
    regs = region_list_to_objects(temp)

    actual = objects_to_ds9_string(regs, coordsys='fk5')
    reference_file = get_pkg_data_filename('data/fk5_reference.reg')
    with open(reference_file, 'r') as fh:
        desired = fh.read()

    assert actual == desired
