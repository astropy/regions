# This file sets up detailed tests for computing masks with reference images

import itertools
import pytest

from regions.core import PixCoord
from regions.shapes.circle import CirclePixelRegion
from regions.shapes.ellipse import EllipsePixelRegion
from regions.shapes.rectangle import RectanglePixelRegion
from astropy.io import fits
from astropy import units as u


REGIONS = [CirclePixelRegion(PixCoord(3.981987, 4.131378), 3.3411),
           EllipsePixelRegion(PixCoord(3.981987, 4.131378), 3.3411, 2.2233, angle=32 * u.deg),
           RectanglePixelRegion(PixCoord(3.981987, 4.131378), 4.3411, 5.2233, angle=32 * u.deg)]

MODES = [{'mode': 'center'},
         {'mode': 'exact'},
         {'mode': 'subpixels', 'subpixels': 1},
         {'mode': 'subpixels', 'subpixels': 5},
         {'mode': 'subpixels', 'subpixels': 18}]


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
