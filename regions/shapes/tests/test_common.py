# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Base class for all shape tests

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_equal, assert_allclose
import pytest

from ...core import PixCoord, BoundingBox
from .utils import ASTROPY_LT_13, HAS_SHAPELY  # noqa


class BaseTestRegion(object):

    def test_repr(self):
        assert repr(self.reg).replace(" ","") == self.expected_repr.replace(" ","")

    def test_str(self):
        assert str(self.reg).replace(" ","") == self.expected_str.replace(" ","")


class BaseTestPixelRegion(BaseTestRegion):

    def test_area(self):
        try:
            assert_allclose(self.reg.area, self.expected_area)
        except ImportError:
            pytest.skip()

    def test_mask_area(self):
        try:
            mask = self.reg.to_mask(mode='exact')
            assert_allclose(np.sum(mask.data), self.expected_area)
        except NotImplementedError:
            try:
                mask = self.reg.to_mask(mode='subpixels', subpixels=30)
                assert_allclose(np.sum(mask.data), self.expected_area, rtol=0.005)
            except NotImplementedError:
                pytest.skip()

    @pytest.mark.skipif('not HAS_SHAPELY')
    def test_contains_compared_to_shapely(self):
        from shapely.geometry import Point
        np.random.seed(12345)
        x = np.random.uniform(self.sample_box[0], self.sample_box[1], 1000)
        y = np.random.uniform(self.sample_box[2], self.sample_box[3], 1000)
        inside = self.reg.contains(PixCoord(x, y))
        reg_shapely = self.reg.to_shapely()
        inside_shapely = [reg_shapely.contains(Point(x[i], y[i])) for i in range(len(x))]
        assert_equal(inside, inside_shapely)

    @pytest.mark.skipif('not HAS_SHAPELY')
    def test_bbox_compared_to_shapely(self):
        reg_shapely = self.reg.to_shapely()
        xmin, ymin, xmax, ymax = reg_shapely.bounds
        bbox_shapely = BoundingBox.from_float(xmin, xmax, ymin, ymax)
        assert self.reg.bounding_box == bbox_shapely

    def test_contains_scalar(self):

        if len(self.inside) > 0:
            pixcoord = PixCoord(*self.inside[0])
            assert self.reg.contains(pixcoord)
            assert pixcoord in self.reg

        if len(self.outside) > 0:
            pixcoord = PixCoord(*self.outside[0])
            assert not self.reg.contains(pixcoord)
            assert pixcoord not in self.reg

    def test_contains_array_1d(self):

        pixcoord = PixCoord(*zip(*(self.inside + self.outside)))

        actual = self.reg.contains(pixcoord)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        with pytest.raises(ValueError) as exc:
            pixcoord in self.reg
        assert 'coord must be scalar' in str(exc)

    def test_contains_array_2d(self):

        x, y = zip(*(self.inside + self.outside))
        pixcoord = PixCoord([x] * 3, [y] * 3)

        actual = self.reg.contains(pixcoord)
        assert actual.shape == (3, len(x))
        assert_equal(actual[:, :len(self.inside)], True)
        assert_equal(actual[:, len(self.inside):], False)


class BaseTestSkyRegion(BaseTestRegion):
    # TODO: here we should add inside/outside tests as above
    pass
