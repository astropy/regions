# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import math

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
from astropy.wcs.utils import pixel_to_skycoord

from ..core import PixCoord, PixelRegion, SkyRegion, Mask, BoundingBox
from .._geometry import elliptical_overlap_grid
from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle


__all__ = ['EllipsePixelRegion', 'EllipseSkyRegion']


class EllipsePixelRegion(PixelRegion):
    """
    An ellipse in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        Center position
    major : float
        Major radius
    minor : float
        Minor radius
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse.
        If set to zero (the default), the major
        axis is lined up with the x axis.

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Ellipse2D
        from astropy.coordinates import Angle
        from regions import PixCoord, EllipsePixelRegion
        import matplotlib.pyplot as plt
        x0, y0 = 15, 10
        a, b = 8, 5
        theta = Angle(30, 'deg')
        e = Ellipse2D(amplitude=100., x_0=x0, y_0=y0, a=a, b=b, theta=theta.radian)
        y, x = np.mgrid[0:20, 0:30]
        fig, ax = plt.subplots(1, 1)
        ax.imshow(e(x, y), origin='lower', interpolation='none', cmap='Greys_r')

        center = PixCoord(x=x0, y=y0)
        reg = EllipsePixelRegion(center=center, major=a, minor=b, angle=theta)
        patch = reg.as_patch(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.show()
    """

    def __init__(self, center, major, minor, angle=0. * u.deg, meta=None,
                 visual=None):
        # TODO: use quantity_input to check that angle is an angle
        self.center = center
        self.major = major
        self.minor = minor
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('major', self.major), ('minor', self.minor),
                             ('angle', self.angle)]

    @property
    def area(self):
        """Region area (float)"""
        return math.pi * self.major * self.minor

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        dx = pixcoord.x - self.center.x
        dy = pixcoord.y - self.center.y
        return (((cos_angle * dx + sin_angle * dy) / self.major) ** 2 +
                ((sin_angle * dx + cos_angle * dy) / self.minor) ** 2 <= 1.)

    def to_shapely(self):
        from shapely import affinity
        ellipse = self.center.to_shapely().buffer(self.minor)
        ellipse = affinity.scale(ellipse, xfact=self.major / self.minor, yfact=1)
        return affinity.rotate(ellipse, self.angle.to(u.deg).value)

    def to_sky(self, wcs):
        # TODO: write a pixel_to_skycoord_scale_angle
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, north_angle = skycoord_to_pixel_scale_angle(center, wcs)
        minor = Angle(self.minor / scale, 'deg')
        major = Angle(self.major / scale, 'deg')
        return EllipseSkyRegion(center, major, minor,
                                angle=self.angle - (north_angle - 90 * u.deg),
                                meta=self.meta, visual=self.visual)

    @property
    def bounding_box(self):
        """
        The minimal bounding box (`~regions.BoundingBox`) enclosing the
        exact elliptical region.
        """

        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        ax = self.major * cos_angle
        ay = self.major * sin_angle
        bx = self.minor * -sin_angle
        by = self.minor * cos_angle
        dx = np.sqrt(ax * ax + bx * bx)
        dy = np.sqrt(ay * ay + by * by)

        xmin = self.center.x - dx
        xmax = self.center.x + dx
        ymin = self.center.y - dy
        ymax = self.center.y + dy

        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):

        # NOTE: assumes this class represents a single circle

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that ellipse is at origin
        xmin = float(bbox.ixmin) - 0.5 - self.center.x
        xmax = float(bbox.ixmax) - 0.5 - self.center.x
        ymin = float(bbox.iymin) - 0.5 - self.center.y
        ymax = float(bbox.iymax) - 0.5 - self.center.y

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = elliptical_overlap_grid(
            xmin, xmax, ymin, ymax, nx, ny,
            self.major, self.minor,
            self.angle.to(u.rad).value,
            use_exact, subpixels,
        )

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Ellipse`).
        """
        from matplotlib.patches import Ellipse
        xy = self.center.x, self.center.y
        width = 2 * self.major
        height = 2 * self.minor
        # From the docstring: MPL expects "rotation in degrees (anti-clockwise)"
        angle = self.angle.to('deg').value
        return Ellipse(xy=xy, width=width, height=height, angle=angle, **kwargs)


class EllipseSkyRegion(SkyRegion):
    """
    An ellipse defined using sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Center position
    major : `~astropy.units.Quantity`
        Major radius
    minor : `~astropy.units.Quantity`
        Minor radius
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse.
        If set to zero (the default), the major
        axis is lined up with the longitude axis of the celestial coordinates.
    """

    def __init__(self, center, major, minor, angle=0. * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that height, width, and angle are angles
        self.center = center
        self.major = major
        self.minor = minor
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('major', self.major), ('minor', self.minor),
                             ('angle', self.angle)]

    def to_pixel(self, wcs):
        center, scale, north_angle = skycoord_to_pixel_scale_angle(self.center, wcs)
        # FIXME: The following line is needed to get a scalar PixCoord
        center = PixCoord(float(center.x), float(center.y))
        minor = self.minor.to('deg').value * scale
        major = self.major.to('deg').value * scale
        return EllipsePixelRegion(center, major, minor,
                                  angle=self.angle + (north_angle - 90 * u.deg),
                                  meta=self.meta, visual=self.visual)
