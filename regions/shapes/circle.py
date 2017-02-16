# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import math

import numpy as np

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs.utils import pixel_to_skycoord

from .polygon import PolygonSkyRegion

from ..core import PixCoord, PixelRegion, SkyRegion, Mask, BoundingBox
from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle
from .._geometry import circular_overlap_grid, rotate_polygon

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        Center position
    radius : float
        Radius
    """

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = PixCoord._validate(center, name='center', expected='scalar')
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('radius', self.radius)]

    @property
    def area(self):
        """Region area (float)."""
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        return self.center.separation(pixcoord) < self.radius

    def to_shapely(self):
        return self.center.to_shapely().buffer(self.radius)

    def to_sky(self, wcs, mode='local', tolerance=None):
        if mode != 'local':
            raise NotImplementedError
        if tolerance is not None:
            raise NotImplementedError

        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        # TODO: this is just called to compute `scale`
        # This is inefficient ... we should have that as a separate function.
        _, scale, _ = skycoord_to_pixel_scale_angle(center, wcs)

        radius = Angle(self.radius / scale, 'deg')
        return CircleSkyRegion(center, radius)

    @property
    def bounding_box(self):
        """
        Bounding box (`~regions.BoundingBox`).
        """
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

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Circle`).
        """
        from matplotlib.patches import Circle
        xy = self.center.x, self.center.y
        radius = self.radius
        return Circle(xy=xy, radius=radius, **kwargs)


class CircleSkyRegion(SkyRegion):
    """
    A circle in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Center position
    radius : `~astropy.units.Quantity`
        Radius in angular units
    """

    def __init__(self, center, radius, meta=None, visual=None):
        # TODO: test that center is a 0D SkyCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('radius', self.radius)]

    @property
    def area(self):
        """Region area (`~astropy.units.Quantity`)"""
        return math.pi * self.radius ** 2

    def contains(self, skycoord):
        return self.center.separation(skycoord) < self.radius

    def to_polygon(self, points=100):
        """
        Convert the circle to a polygon.

        Parameters
        ----------
        points : int, optional
            Number of points in the final polygon.
        """

        # TODO: avoid converting to unit spherical or spherical if already
        #       using a spherical representation

        # Extract longitude/latitude, either from a tuple of two quantities, or
        # a single 2-element Quantity.
        rep = self.center.represent_as('unitspherical')
        longitude, latitude = rep.lon, rep.lat

        # Start off by generating the circle around the North pole
        lon = np.linspace(0., 2 * np.pi, points + 1)[:-1] * u.radian
        lat = np.repeat(0.5 * np.pi - self.radius.to(u.radian).value, points) * u.radian

        # Now rotate it to the correct longitude/latitude
        lon, lat = rotate_polygon(lon, lat, longitude, latitude)

        # Make a new SkyCoord
        vertices_sky = SkyCoord(lon, lat, frame=self.center)

        return PolygonSkyRegion(vertices_sky)

    def to_pixel(self, wcs, mode='local', tolerance=100):
        """
        Given a WCS, convert the circle to a best-approximation circle in pixel
        dimensions.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            A world coordinate system
        mode : 'local' or not
            not implemented
        tolerance : int
            not implemented

        Returns
        -------
        CirclePixelRegion
        """

        if mode == 'local':
            center, scale, _ = skycoord_to_pixel_scale_angle(self.center, wcs)
            # The following line is needed to get a scalar PixCoord
            center = PixCoord(float(center.x), float(center.y))
            radius = self.radius.to('deg').value * scale
            return CirclePixelRegion(center, radius)

        elif mode == 'affine':
            raise NotImplementedError()

        elif mode == 'full':

            return self.to_polygon(points=tolerance).to_pixel(wcs)

        else:
            raise ValueError('mode should be one of local/affine/full')
