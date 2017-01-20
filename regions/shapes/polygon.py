# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from astropy.wcs.utils import skycoord_to_pixel

from ..core import PixelRegion, SkyRegion, Mask, BoundingBox, PixCoord
from .._geometry import polygonal_overlap_grid
from .._geometry.pnpoly import points_in_polygon

__all__ = ['PolygonPixelRegion', 'PolygonSkyRegion']


class PolygonPixelRegion(PixelRegion):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    vertices : `~regions.PixCoord`
        The vertices of the polygon
    """

    def __init__(self, vertices, meta=None, visual=None):
        # TODO: test that vertices is a 1D PixCoord
        self.vertices = vertices
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('vertices', self.vertices)]

    @property
    def area(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            x = np.array([pixcoord.x], dtype=float)
            y = np.array([pixcoord.y], dtype=float)
        else:
            x = pixcoord.x.astype(float)
            y = pixcoord.y.astype(float)
        return points_in_polygon(x, y,
                                 self.vertices.x.astype(float),
                                 self.vertices.y.astype(float)).astype(bool)

    def to_shapely(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_sky(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError

    @property
    def bounding_box(self, mode='center'):
        xmin = self.vertices.x.min()
        xmax = self.vertices.x.max()
        ymin = self.vertices.y.min()
        ymax = self.vertices.y.max()
        return BoundingBox._from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5
        xmax = float(bbox.ixmax) - 0.5
        ymin = float(bbox.iymin) - 0.5
        ymax = float(bbox.iymax) - 0.5

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = polygonal_overlap_grid(
            xmin, xmax, ymin, ymax, nx, ny,
            self.vertices.x.astype(float),
            self.vertices.y.astype(float),
            use_exact, subpixels)

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """
        Matplotlib patch object for this region (`matplotlib.patches.Polygon`).
        """
        from matplotlib.patches import Polygon
        xy = np.vstack([self.vertices.x, self.vertices.y]).transpose()
        return Polygon(xy=xy, **kwargs)


class PolygonSkyRegion(SkyRegion):
    """
    A polygon in sky coordinates.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon
    """

    def __init__(self, vertices, meta=None, visual=None):
        # TODO: test that vertices is a 1D SkyCoord
        self.vertices = vertices
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('vertices', self.vertices)]

    @property
    def area(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def contains(self, skycoord):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_pixel(self, wcs, mode='local', tolerance=None):

        if mode == 'local':
            x, y = skycoord_to_pixel(self.vertices, wcs)
            vertices_pix = PixCoord(x, y)
            return PolygonPixelRegion(vertices_pix)
        else:

            raise NotImplementedError()
