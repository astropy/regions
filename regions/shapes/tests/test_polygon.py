# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy.coordinates import SkyCoord

from ..._utils.examples import make_example_dataset
from ...core import PixCoord, BoundingBox
from ...tests.helpers import make_simple_wcs
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from .utils import HAS_MATPLOTLIB  # noqa
from .test_common import BaseTestPixelRegion, BaseTestSkyRegion


@pytest.fixture(scope='session')
def wcs():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs


class TestPolygonPixelRegion(BaseTestPixelRegion):

    reg = PolygonPixelRegion(PixCoord([1, 3, 1], [1, 1, 4]))
    sample_box = [0, 4, 0, 5]
    inside = [(2, 2)]
    outside = [(3, 2), (3, 3)]
    expected_area = 3
    expected_repr = '<PolygonPixelRegion(vertices=PixCoord(x=[1 3 1], y=[1 1 4]))>'
    expected_str = ('Region: PolygonPixelRegion\nvertices: PixCoord(x=[1 3 1],'
                    ' y=[1 1 4])')

    # We will be using this polygon for basic tests:
    #
    # (3,0) *
    #       |\
    #       | \
    #       |  \
    # (0,0) *---* (2, 0)
    #

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.vertices.x, [1, 3, 1])
        assert_allclose(reg.vertices.y, [1, 1, 4])
        assert reg.visual == {}
        assert reg.meta == {}

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.vertices.x, self.reg.vertices.x)
        assert_allclose(reg_new.vertices.y, self.reg.vertices.y)

    def test_bounding_box(self):
        bbox = self.reg.bounding_box
        assert bbox == BoundingBox(ixmin=1, ixmax=4, iymin=1, iymax=5)

    def test_to_mask(self):
        # The true area of this polygon is 3
        # We only do very low-precision asserts below,
        # because results can be unstable with points
        # on the edge of the polygon.
        # Basically we check that mask.data is filled
        # with something useful at all.

        # Bounding box and output shape is independent of subpixels,
        # so we only assert on it once here, not in the other cases below
        mask = self.reg.to_mask(mode='center', subpixels=1)
        assert 2 <= np.sum(mask.data) <= 6
        assert mask.bbox == BoundingBox(ixmin=1, ixmax=4, iymin=1, iymax=5)
        assert mask.data.shape == (4, 3)

        # Test more cases for to_mask
        # This example is with the default: subpixels=5
        mask = self.reg.to_mask(mode='subpixels')
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=8)
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=9)
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=10)
        assert 2 <= np.sum(mask.data) <= 6

        with pytest.raises(NotImplementedError):
            self.reg.to_mask(mode='exact')

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        expected = [[1, 1], [3, 1], [1, 4], [1, 1]]
        assert_allclose(patch.xy, expected)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(3, 1), -90 * u.deg)
        assert_allclose(reg.vertices.x, [3, 3, 6])
        assert_allclose(reg.vertices.y, [3, 1, 3])


class TestPolygonSkyRegion(BaseTestSkyRegion):

    reg = PolygonSkyRegion(SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg))

    expected_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                'dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,  '
                '4.)]>)>')
    expected_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
               ' (ra, dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,'
               '  4.)]>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.vertices.ra.deg, [3, 4, 3])
        assert reg.visual == {}
        assert reg.meta == {}

    def test_transformation(self, wcs):

        pixpoly = self.reg.to_pixel(wcs)

        assert_allclose(pixpoly.vertices.x, [11.187992, 10.976332, 11.024032], atol=1e-5)
        assert_allclose(pixpoly.vertices.y, [1.999486, 2.039001, 2.077076], atol=1e-5)

        poly = pixpoly.to_sky(wcs)

        # TODO: we should probably assert something about frame attributes,
        # or generally some better way to check if two SkyCoord are the same?
        # For now, we use the folloing line to transform back to ICRS (`poly` is in Galactic, same as WCS)
        poly = PolygonSkyRegion(vertices=poly.vertices.transform_to(self.reg.vertices))
        assert_quantity_allclose(poly.vertices.data.lon, self.reg.vertices.data.lon, atol=1e-3 * u.deg)
        assert_quantity_allclose(poly.vertices.data.lat, self.reg.vertices.data.lat, atol=1e-3 * u.deg)

    def test_contains(self, wcs):
        position = SkyCoord([1, 3.25] * u.deg, [2, 3.75] * u.deg)
        # 1,2 is outside, 3.25,3.75 should be inside the triangle...
        assert all(self.reg.contains(position, wcs) == np.array([False, True], dtype='bool'))
