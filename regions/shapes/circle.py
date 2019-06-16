# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import math

from astropy.coordinates import Angle
from astropy.wcs.utils import pixel_to_skycoord

from ..core import PixCoord, PixelRegion, SkyRegion, RegionMask, BoundingBox
from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle
from .._geometry import circular_overlap_grid
from ..core.attributes import (ScalarSky, ScalarPix, QuantityLength,
                               ScalarLength, RegionVisual, RegionMeta)

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle defined using pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        Center position
    radius : `float`
        Radius
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

    Examples
    --------

    .. plot::
        :include-source:

        from regions import PixCoord, CirclePixelRegion
        import matplotlib.pyplot as plt

        x, y = 6, 6
        radius = 5.5

        fig, ax = plt.subplots(1, 1)

        center = PixCoord(x=x, y=y)
        reg = CirclePixelRegion(center=center, radius=radius)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(0, 15)
        plt.ylim(0, 15)
        ax.set_aspect('equal')
    """

    center = ScalarPix('center')
    radius = ScalarLength('radius')

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('radius',)

    @property
    def area(self):
        """Region area (`float`)."""
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        in_circle = self.center.separation(pixcoord) < self.radius
        if self.meta.get('include', True):
            return in_circle
        else:
            return np.logical_not(in_circle)

    def to_sky(self, wcs):
        # TODO: write a pixel_to_skycoord_scale_angle
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, _ = skycoord_to_pixel_scale_angle(center, wcs)
        radius = Angle(self.radius / scale, 'deg')
        return CircleSkyRegion(center, radius, self.meta, self.visual)

    @property
    def bounding_box(self):
        """Bounding box (`~regions.BoundingBox`)."""
        xmin = self.center.x - self.radius
        xmax = self.center.x + self.radius
        ymin = self.center.y - self.radius
        ymax = self.center.y + self.radius

        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=1):
        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5 - self.center.x
        xmax = float(bbox.ixmax) - 0.5 - self.center.x
        ymin = float(bbox.iymin) - 0.5 - self.center.y
        ymax = float(bbox.iymax) - 0.5 - self.center.y

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = circular_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                         self.radius, use_exact, subpixels)

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Matplotlib patch object for this region (`matplotlib.patches.Circle`)

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.patches.Circle` object accepts

        Returns
        -------
        patch : `~matplotlib.patches.Circle`
            Matplotlib circle patch
        """
        from matplotlib.patches import Circle
        xy = self.center.x - origin[0], self.center.y - origin[1]
        radius = self.radius

        mpl_params = self.mpl_properties_default('patch')
        mpl_params.update(kwargs)

        return Circle(xy=xy, radius=radius, **mpl_params)

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
        region : `CirclePixelRegion`
            Rotated region (an independent copy)
        """
        center = self.center.rotate(center, angle)
        return self.copy(center=center)


class CircleSkyRegion(SkyRegion):
    """
    A circle defined using sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Center position
    radius : `~astropy.units.Quantity`
        Radius in angular units
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    center = ScalarSky('center')
    radius = QuantityLength("radius")

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = ('radius',)

    def to_pixel(self, wcs):
        center, scale, _ = skycoord_to_pixel_scale_angle(self.center, wcs)
        radius = self.radius.to('deg').value * scale
        return CirclePixelRegion(center, radius, self.meta, self.visual)
