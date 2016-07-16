# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_allclose
from ..examples import make_example_wcs


def test_make_example_wcs():
    wcs = make_example_wcs()
    assert_allclose(wcs.wcs.crval, (0, 0))
    assert_allclose(wcs.wcs.crpix, (18, 9))
    assert_allclose(wcs.wcs.cdelt, (10, 10))
    assert wcs.wcs.ctype[0] == 'GLON-AIT'
    assert wcs.wcs.ctype[1] == 'GLAT-AIT'
