# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file sets up detailed tests for computing masks with reference images.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import itertools
import pytest

from astropy import units as u

from ...core import PixCoord
from ...shapes.circle import CirclePixelRegion
from ...shapes.ellipse import EllipsePixelRegion
from ...shapes.rectangle import RectanglePixelRegion
from ...shapes.polygon import PolygonPixelRegion

REGIONS = [
    CirclePixelRegion(PixCoord(3.981987, 4.131378), radius=3.3411),
    EllipsePixelRegion(PixCoord(5.981987, 4.131378), width=10.4466, height=6.6822, angle=12 * u.deg),
    RectanglePixelRegion(PixCoord(5.981987, 4.131378), width=7.2233, height=4.3411, angle=12 * u.deg),
    PolygonPixelRegion(PixCoord([-2.334, 3.631, 1.122, -1.341], [-3.121, -2.118, 2.987, 1.759])),
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
    elif isinstance(value, PolygonPixelRegion):
        return 'poly'
    else:
        return '-'.join('{0}_{1}'.format(key, value) for key, value in sorted(value.items()))

# There's a bug in numpy: `numpy.savetxt` doesn't accept unicode
# on Python 2. Bytes works on Python 2 and 3, so we're using that here.
# https://github.com/numpy/numpy/pull/4053#issuecomment-263808221
@pytest.mark.array_compare(
    fmt='text',
    write_kwargs={'fmt': b'%12.8e'},
)
@pytest.mark.parametrize(
    ('region', 'mode'),
    itertools.product(REGIONS, MODES),
    ids=label,
)
def test_to_mask(region, mode):
    try:
        mask = region.to_mask(**mode)
    except NotImplementedError:
        pytest.xfail()
    return mask.data.astype(float)
