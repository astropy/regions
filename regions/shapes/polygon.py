# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np

from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

from ..core import PixelRegion, SkyRegion, RegionMask, BoundingBox, PixCoord
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
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        from astropy.coordinates import Angle
        from regions import PixCoord, PolygonPixelRegion
        import matplotlib.pyplot as plt

        x, y = [45, 45, 55, 60], [75, 70, 65, 75]
        fig, ax = plt.subplots(1, 1)

        vertices = PixCoord(x=x, y=y)
        reg = PolygonPixelRegion(vertices=vertices)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(30, 80)
        plt.ylim(50, 80)
        ax.set_aspect('equal')
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
        raise NotImplementedError

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, 'pixcoord')
        x = np.atleast_1d(np.asarray(pixcoord.x, dtype=float))
        y = np.atleast_1d(np.asarray(pixcoord.y, dtype=float))
        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        shape = x.shape
        mask = points_in_polygon(x.flatten(), y.flatten(), vx, vy).astype(bool)
        in_poly = mask.reshape(shape)
        if self.meta.get('include', True):
            return in_poly
        else:
            return np.logical_not(in_poly)

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

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Matplotlib patch object for this region (`matplotlib.patches.Polygon`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.patches.Polygon` object accepts

        Returns
        -------
        patch : `~matplotlib.patches.Polygon`
            Matplotlib polygon patch
        """
        from matplotlib.patches import Polygon
        xy = np.vstack([self.vertices.x - origin[0],
                        self.vertices.y - origin[1]]).transpose()

        mpl_params = self.mpl_properties_default('patch')
        mpl_params.update(kwargs)

        return Polygon(xy=xy, **mpl_params)

    def rotate(self, center, angle):
        """Make a rotated region.

        Rotates counter-clockwise for positive ``angle``.

        Parameters
        ----------
        center : `PixCoord`
            Rotation center point
        angle : `~astropy.coordinates.Angle`
            Rotation angle

        Returns
        -------
        region : `PolygonPixelRegion`
            Rotated region (an independent copy)
        """
        vertices = self.vertices.rotate(center, angle)
        return self.copy(vertices=vertices)


class PolygonSkyRegion(SkyRegion):
    """
    A polygon defined using vertices in sky coordinates.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
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
