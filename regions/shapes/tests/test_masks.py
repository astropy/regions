# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file sets up detailed tests for computing masks with reference images.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import itertools
from astropy.tests.helper import pytest
from astropy import units as u
from ...core import PixCoord
from ...shapes.circle import CirclePixelRegion
from ...shapes.ellipse import EllipsePixelRegion
from ...shapes.rectangle import RectanglePixelRegion

REGIONS = [
    CirclePixelRegion(PixCoord(3.981987, 4.131378), radius=3.3411),
    EllipsePixelRegion(PixCoord(3.981987, 4.131378), major=2.2233, minor=3.3411, angle=32 * u.deg),
    RectanglePixelRegion(PixCoord(3.981987, 4.131378), width=5.2233, height=4.3411, angle=32 * u.deg),
]

MODES = [
    {'mode': 'center'},
    {'mode': 'exact'},
    {'mode': 'subpixels', 'subpixels': 1},
    {'mode': 'subpixels', 'subpixels': 5},
    {'mode': 'subpixels', 'subpixels': 18},
]


def label(value):
    if isinstance(value, CirclePixelRegion):
        return 'circ'
    elif isinstance(value, EllipsePixelRegion):
        return 'elli'
    elif isinstance(value, RectanglePixelRegion):
        return 'rect'
    else:
        return '-'.join('{0}_{1}'.format(key, value) for key, value in sorted(value.items()))


@pytest.mark.array_compare(fmt='text', write_kwargs={'fmt': '%12.8e'})
@pytest.mark.parametrize(('region', 'mode'), itertools.product(REGIONS, MODES), ids=label)
def test_to_mask(region, mode):
    try:
        mask = region.to_mask(**mode)
    except NotImplementedError:
        pytest.xfail()
    return mask.data.astype(float)
