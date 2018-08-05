# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

from ..core import PixelRegion, SkyRegion, Mask, BoundingBox, PixCoord
from .._geometry import polygonal_overlap_grid
from .._geometry.pnpoly import points_in_polygon
from ..core.attributes import OneDPix, OneDSky, RegionMeta, RegionVisual

__all__ = ['PolygonPixelRegion', 'PolygonSkyRegion']


class PolygonPixelRegion(PixelRegion):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    vertices : `~regions.PixCoord`
        The vertices of the polygon
    meta: `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual: `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    vertices = OneDPix('vertices')

    def __init__(self, vertices, meta=None, visual=None):
        self.vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('vertices',)

    @property
    def area(self):
        """Region area (float)."""
        # FIXME: for now we use shapely, but could try and avoid the dependency in future
        return self.to_shapely().area

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, 'pixcoord')
        x = np.atleast_1d(np.asarray(pixcoord.x, dtype=float))
        y = np.atleast_1d(np.asarray(pixcoord.y, dtype=float))
        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        shape = x.shape
        mask = points_in_polygon(x.flatten(), y.flatten(), vx, vy).astype(bool)
        return mask.reshape(shape)

    def to_shapely(self):
        from shapely.geometry import Polygon
        return Polygon(list(zip(self.vertices.x, self.vertices.y)))

    def to_sky(self, wcs):
        vertices_sky = pixel_to_skycoord(self.vertices.x, self.vertices.y, wcs)
        return PolygonSkyRegion(vertices=vertices_sky)

    @property
    def bounding_box(self):
        xmin = self.vertices.x.min()
        xmax = self.vertices.x.max()
        ymin = self.vertices.y.min()
        ymax = self.vertices.y.max()
        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5
        xmax = float(bbox.ixmax) - 0.5
        ymin = float(bbox.iymin) - 0.5
        ymax = float(bbox.iymax) - 0.5

        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        fraction = polygonal_overlap_grid(
            xmin, xmax, ymin, ymax,
            nx, ny, vx, vy, use_exact, subpixels,
        )

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Polygon`)."""
        from matplotlib.patches import Polygon
        xy = np.vstack([self.vertices.x, self.vertices.y]).transpose()
        return Polygon(xy=xy, **kwargs)


class PolygonSkyRegion(SkyRegion):
    """
    A polygon defined using vertices in sky coordinates.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon
    meta: `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual: `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    vertices = OneDSky('vertices')

    def __init__(self, vertices, meta=None, visual=None):
        self.vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('vertices',)

    def to_pixel(self, wcs):
        x, y = skycoord_to_pixel(self.vertices, wcs)
        vertices_pix = PixCoord(x, y)
        return PolygonPixelRegion(vertices_pix)
