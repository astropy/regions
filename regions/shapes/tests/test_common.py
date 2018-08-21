# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Base class for all shape tests

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from numpy.testing import assert_equal, assert_allclose
import pytest

from ..polygon import PolygonPixelRegion
from ...core import PixCoord, BoundingBox

class BaseTestRegion(object):

    def test_repr(self):
        assert repr(self.reg).replace(" ","") == self.expected_repr.replace(" ","")

    def test_str(self):
        assert str(self.reg).replace(" ","") == self.expected_str.replace(" ","")


class BaseTestPixelRegion(BaseTestRegion):

    def test_area(self):
        # TODO: remove the pytest.skip once polygon area is implemented
        if isinstance(self.reg, PolygonPixelRegion):
            pytest.skip()

        assert_allclose(self.reg.area, self.expected_area)

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
