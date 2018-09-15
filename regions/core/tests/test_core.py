# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
from ...shapes import CircleSkyRegion, PointSkyRegion

try:
    import matplotlib

    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_sky_region_plot():
    assert 'TODO'
