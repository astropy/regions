from __future__ import absolute_import, division, print_function, unicode_literals
from ..ds9_language import ds9_parser
from astropy.utils.data import get_pkg_data_filename

def test_fk5():
    filename = get_pkg_data_filename('data/ds9.fk5.reg')
    regs = ds9_parser(filename)
    
