# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from astropy.tests.helper import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from ...core import PixCoord, BoundingBox
from ..polygon import PolygonPixelRegion, PolygonSkyRegion
from .utils import ASTROPY_LT_13, HAS_MATPLOTLIB


class TestPolygonPixelRegion:
    def setup(self):
        # We will be using this polygon for basic tests:
        #
        # (3,0) *
        #       |\
        #       | \
        #       |  \
        # (0,0) *---* (2, 0)
        #
        vertices = PixCoord([0, 2, 0], [0, 0, 3])
        self.reg = PolygonPixelRegion(vertices)
        self.pixcoord_inside = PixCoord(1, 1)
        self.pixcoord_outside = PixCoord(1, 2)

    def test_repr_str(self):
        reg_repr = ('<PolygonPixelRegion(vertices=PixCoord(x=[0 2 0], '
                    'y=[0 0 3]))>')
        assert repr(self.reg) == reg_repr

        reg_str = ('Region: PolygonPixelRegion\nvertices: PixCoord(x=[0 2 0],'
                   ' y=[0 0 3])')
        assert str(self.reg) == reg_str

    def test_contains_scalar(self):
        assert self.reg.contains(self.pixcoord_inside)
        assert self.pixcoord_inside in self.reg

        assert not self.reg.contains(self.pixcoord_outside)
        assert self.pixcoord_outside not in self.reg

    def test_contains_array_1d(self):
        pixcoord = PixCoord([1, 1], [1, 2])
        actual = self.reg.contains(pixcoord)
        expected = [True, False]
        assert_equal(actual, expected)

        with pytest.raises(ValueError) as exc:
            pixcoord in self.reg
        assert 'coord must be scalar' in str(exc)

    def test_contains_array_2d(self):
        pixcoord = PixCoord(
            [[1, 1, 1], [1, 1, 1]],
            [[1, 1, 1], [2, 2, 2]],
        )
        actual = self.reg.contains(pixcoord)
        expected = [[True, True, True], [False, False, False]]
        assert_equal(actual, expected)

    def test_bounding_box(self):
        bbox = self.reg.bounding_box
        assert bbox == BoundingBox(ixmin=0, ixmax=3, iymin=0, iymax=4)

    def test_to_mask(self):
        # The true area of this polygon is 3

        mask = self.reg.to_mask(mode='center', subpixels=1)
        assert_allclose(np.sum(mask.data), 5)

        # Bounding box and output shape is independent of subpixels,
        # so we only assert on it once here, not in the other cases below
        assert mask.bbox == BoundingBox(ixmin=0, ixmax=3, iymin=0, iymax=4)
        assert mask.data.shape == (4, 3)

        # This example is with the default: subpixels=5
        mask = self.reg.to_mask(mode='subpixels')
        assert_allclose(np.sum(mask.data), 3.48)

        mask = self.reg.to_mask(mode='subpixels', subpixels=8)
        assert_allclose(np.sum(mask.data), 3.0)

        mask = self.reg.to_mask(mode='subpixels', subpixels=9)
        assert_allclose(np.sum(mask.data), 2.6790123456790127)

        mask = self.reg.to_mask(mode='subpixels', subpixels=10)
        assert_allclose(np.sum(mask.data), 3.0)

        with pytest.raises(NotImplementedError):
            self.reg.to_mask(mode='exact')

    @pytest.mark.skipif('not HAS_MATPLOTLIB')
    def test_as_patch(self):
        patch = self.reg.as_patch()
        expected = [[0, 0], [2, 0], [0, 3], [0, 0]]
        assert_allclose(patch.xy, expected)


class TestPolygonSkyRegion:
    def setup(self):
        vertices = SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg)
        self.poly = PolygonSkyRegion(vertices)

    def test_repr_str(self):
        if ASTROPY_LT_13:
            reg_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                        'dec) in deg\n    [(3.0, 3.0), (4.0, 4.0), (3.0, '
                        '4.0)]>)>')
            reg_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
                       ' (ra, dec) in deg\n    [(3.0, 3.0), (4.0, 4.0), (3.0,'
                       ' 4.0)]>')
        else:
            reg_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                        'dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,  '
                        '4.)]>)>')
            reg_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
                       ' (ra, dec) in deg\n    [( 3.,  3.), ( 4.,  4.), ( 3.,'
                       '  4.)]>')

        assert repr(self.poly) == reg_repr
        assert str(self.poly) == reg_str
