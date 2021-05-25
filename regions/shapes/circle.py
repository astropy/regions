# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines circular regions in both pixel and sky coordinates.
"""

import math

from astropy.coordinates import Angle
import astropy.units as u
from astropy.wcs.utils import pixel_to_skycoord
import numpy as np

from ..core.attributes import (ScalarPix, ScalarLength, QuantityLength,
                               ScalarSky)
from ..core.bounding_box import BoundingBox
from ..core.core import PixelRegion, SkyRegion
from ..core.mask import RegionMask
from ..core.metadata import RegionMeta, RegionVisual
from ..core.pixcoord import PixCoord
from .._utils.wcs_helpers import pixel_scale_angle_at_skycoord
from .._geometry import circular_overlap_grid

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle defined using pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The center position.
    radius : float
        The radius.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.

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

    _params = ('center', 'radius')
    center = ScalarPix('center')
    radius = ScalarLength('radius')

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
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
        _, pixscale, _ = pixel_scale_angle_at_skycoord(center, wcs)
        radius = Angle(self.radius * u.pix * pixscale, 'arcsec')
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

        # Find position of pixel edges and recenter so that circle is at
        # origin
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
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Circle`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Circle`. These keywords will override
            any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Circle`
            A matplotlib circle patch.
        """
        from matplotlib.patches import Circle

        xy = self.center.x - origin[0], self.center.y - origin[1]
        radius = self.radius

        mpl_kwargs = self._define_mpl_kwargs(artist='Patch')
        mpl_kwargs.update(kwargs)

        return Circle(xy=xy, radius=radius, **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Postive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `CirclePixelRegion`
            The rotated region (which is an independent copy).
        """
        center = self.center.rotate(center, angle)
        return self.copy(center=center)


class CircleSkyRegion(SkyRegion):
    """
    A circle defined using sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The center position.
    radius : `~astropy.units.Quantity`
        The radius in angular units.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('center', 'radius')
    center = ScalarSky('center')
    radius = QuantityLength("radius")

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        center, pixscale, _ = pixel_scale_angle_at_skycoord(self.center, wcs)
        radius = (self.radius / pixscale).to(u.pix).value
        return CirclePixelRegion(center, radius, self.meta, self.visual)
