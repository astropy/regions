# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines circular regions in both pixel and sky coordinates.
"""

import math

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle

from regions._geometry import circular_overlap_grid
from regions._utils.wcs_helpers import pixel_scale_angle_at_skycoord
from regions.core.attributes import (PositiveScalar, PositiveScalarAngle,
                                     RegionMetaDescr, RegionVisualDescr,
                                     ScalarPixCoord, ScalarSkyCoord)
from regions.core.bounding_box import RegionBoundingBox
from regions.core.core import PixelRegion, SkyRegion
from regions.core.mask import RegionMask
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord


__all__ = ['CirclePixelRegion', 'CircleSkyRegion', 'CircleSectorPixelRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle defined using pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The center position.
    radius : float
        The radius in pixels.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, CirclePixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = CirclePixelRegion(PixCoord(x=8, y=7), radius=3.5)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Circle')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_aspect('equal')
    """

    _params = ('center', 'radius')
    _mpl_artist = 'Patch'
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    radius = PositiveScalar('The radius in pixels as a float.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

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
        center = wcs.pixel_to_world(self.center.x, self.center.y)
        _, pixscale, _ = pixel_scale_angle_at_skycoord(center, wcs)
        radius = Angle(self.radius * u.pix * pixscale, 'arcsec')
        return CircleSkyRegion(center, radius, meta=self.meta.copy(),
                               visual=self.visual.copy())

    @property
    def bounding_box(self):
        """
        Bounding box (`~regions.RegionBoundingBox`).
        """
        xmin = self.center.x - self.radius
        xmax = self.center.x + self.radius
        ymin = self.center.y - self.radius
        ymax = self.center.y + self.radius

        return RegionBoundingBox.from_float(xmin, xmax, ymin, ymax)

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

        use_exact = 0 if mode == 'subpixels' else 1

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
        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Circle(xy=xy, radius=radius, **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

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
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'radius')
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    radius = PositiveScalarAngle('The radius as a |Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        center, pixscale, _ = pixel_scale_angle_at_skycoord(self.center, wcs)
        radius = (self.radius / pixscale).to(u.pix).value
        return CirclePixelRegion(center, radius, meta=self.meta.copy(),
                                 visual=self.visual.copy())


class CircleSectorPixelRegion(PixelRegion):
    """
    A circle sector defined using pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The center position.
    radius : float
        The radius in pixels.
    angle_start: `~astropy.units.Quantity`, optional
        The start angle of the sector, measured anti-clockwise.
    angle_stop : `~astropy.units.Quantity`, optional
        The stop angle of the sector, measured anti-clockwise.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from astropy import units as u
        from regions import PixCoord, CircleSectorPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = CircleSectorPixelRegion(PixCoord(x=8, y=7), radius=3.5, angle_start=0 * u.deg,
                                      angle_stop=120 * u.deg)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Circle')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_aspect('equal')
    """

    _params = ('center', 'radius', 'angle_start', 'angle_stop')
    _mpl_artist = 'Patch'
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    radius = PositiveScalar('The radius in pixels as a float.')
    angle_start = ScalarAngle('The start angle measured anti-clockwise as a '
                              '|Quantity| angle.')
    angle_stop = ScalarAngle('The stop angle measured anti-clockwise as a '
                             '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, radius, angle_start=0 * u.deg, angle_stop=360 * u.deg,
                 meta=None, visual=None):
        self.center = center
        self.radius = radius

        if angle_start >= angle_stop:
            raise ValueError('angle_stop must be greater than angle_start')

        self.angle_start = angle_start
        self.angle_stop = angle_stop
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def theta(self):
        """Opening angle of the sector (`~astropy.coordinates.Angle`)"""
        return self.angle_stop - self.angle_start

    @property
    def area(self):
        return self.radius ** 2 * self.theta.to_value("rad") / 2.

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        in_circle = self.center.separation(pixcoord) < self.radius

        dx = pixcoord.x - self.center.x
        dy = pixcoord.y - self.center.y
        angle = (Angle(np.arctan2(dy, dx), "rad") - self.angle_start).wrap_at("360d")

        in_angle = (angle > 0 * u.deg) & (angle < self.theta)
        in_sector = in_circle & in_angle

        if self.meta.get('include', True):
            return in_sector
        else:
            return np.logical_not(in_sector)

    def to_sky(self, wcs):
        raise NotImplementedError

    def to_mask(self, **kwargs):
        raise NotImplementedError

    @property
    def bounding_box(self):
        """Bounding box (`~regions.RegionBoundingBox`)."""
        x_start = self.radius * np.cos(self.angle_start)
        y_start = self.radius * np.sin(self.angle_start)

        x_stop = self.radius * np.cos(self.angle_stop)
        y_stop = self.radius * np.sin(self.angle_stop)

        def wrap(angle):
            return Angle(angle).wrap_at("360d")

        cross_0 = wrap(self.angle_start) > wrap(self.angle_stop)
        cross_90 = wrap(self.angle_start - 90 * u.deg) > wrap(self.angle_stop - 90 * u.deg)
        cross_180 = wrap(self.angle_start - 180 * u.deg) > wrap(self.angle_stop - 180 * u.deg)
        cross_270 = wrap(self.angle_start - 270 * u.deg) > wrap(self.angle_stop - 270 * u.deg)

        xmin = self.center.x + min(np.where(cross_180, -self.radius, min(x_start, x_stop)), 0)
        xmax = self.center.x + max(np.where(cross_0, self.radius, max(x_start, x_stop)), 0)
        ymin = self.center.y + min(np.where(cross_270, -self.radius, min(y_start, y_stop)), 0)
        ymax = self.center.y + max(np.where(cross_90, self.radius, max(y_start, y_stop)), 0)

        return RegionBoundingBox.from_float(xmin, xmax, ymin, ymax)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Wedge).

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
        artist : `~matplotlib.patches.Wedge`
            A matplotlib circle patch.
        """
        from matplotlib.patches import Wedge

        center = self.center.x - origin[0], self.center.y - origin[1]
        radius = self.radius
        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Wedge(center=center, r=radius, theta1=self.angle_start.to_value("deg"),
                     theta2=self.angle_stop.to_value("deg"), **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

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
        angle_start = self.angle_start + angle
        angle_stop = self.angle_stop + angle
        return self.copy(center=center, angle_start=angle_start, angle_stop=angle_stop)
