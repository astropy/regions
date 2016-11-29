# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import math

from astropy.coordinates import Angle
from astropy.wcs.utils import pixel_to_skycoord

from ..core import PixelRegion, SkyRegion, Mask, BoundingBox
from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle
from .._geometry import circular_overlap_grid

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
        # TODO: test that center is a 0D PixCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            radius=self.radius,
        )
        fmt = '{name}\ncenter: {center}\nradius: {radius}'
        return fmt.format(**data)

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
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

        # Find exact bounds
        xmin = self.center.x - self.radius
        xmax = self.center.x + self.radius
        ymin = self.center.y - self.radius
        ymax = self.center.y + self.radius

        return BoundingBox._from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=1):

        # NOTE: assumes this class represents a single circle

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

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            radius=self.radius,
        )
        fmt = '{name}\ncenter: {center}\nradius: {radius}'
        return fmt.format(**data)

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, skycoord):
        return self.center.separation(skycoord) < self.radius

    def to_pixel(self, wcs, mode='local', tolerance=None):
        """
        Given a WCS, convert the circle to a best-approximation circle in pixel
        dimensions.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            A world coordinate system
        mode : 'local' or not
            not implemented
        tolerance : None
            not implemented

        Returns
        -------
        CirclePixelRegion
        """

        if mode != 'local':
            raise NotImplementedError
        if tolerance is not None:
            raise NotImplementedError

        center, scale, _ = skycoord_to_pixel_scale_angle(self.center, wcs)
        radius = self.radius.to('deg').value * scale

        return CirclePixelRegion(center, radius)
