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
        The position of the center of the ellipse.
    width : float
        The width of the ellipse (before rotation) in pixels
    height : float
        The height of the ellipse (before rotation) in pixels
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the x axis.

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
        reg = EllipsePixelRegion(center=center, width=a, height=b, angle=theta)
        patch = reg.as_patch(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.show()
    """

    def __init__(self, center, width, height, angle=0. * u.deg, meta=None,
                 visual=None):
        # TODO: use quantity_input to check that angle is an angle
        self.center = PixCoord._validate(center, name='center', expected='scalar')
        self.width = width
        self.height = height
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('width', self.width), ('height', self.height),
                             ('angle', self.angle)]

    @property
    def area(self):
        """Region area (float)"""
        return math.pi * self.width * self.height * 0.25

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        dx = pixcoord.x - self.center.x
        dy = pixcoord.y - self.center.y
        return ((2 * (cos_angle * dx + sin_angle * dy) / self.width) ** 2 +
                (2 * (sin_angle * dx - cos_angle * dy) / self.height) ** 2 <= 1.)

    def to_shapely(self):
        from shapely import affinity
        ellipse = self.center.to_shapely().buffer(0.5 * self.height)
        ellipse = affinity.scale(ellipse, xfact=self.width / self.height, yfact=1)
        return affinity.rotate(ellipse, self.angle.to(u.deg).value)

    def to_sky(self, wcs):
        # TODO: write a pixel_to_skycoord_scale_angle
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, north_angle = skycoord_to_pixel_scale_angle(center, wcs)
        height = Angle(self.height / scale, 'deg')
        width = Angle(self.width / scale, 'deg')
        return EllipseSkyRegion(center, width, height,
                                angle=self.angle - (north_angle - 90 * u.deg),
                                meta=self.meta, visual=self.visual)

    @property
    def bounding_box(self):
        """
        The minimal bounding box (`~regions.BoundingBox`) enclosing the
        exact elliptical region.
        """

        # We use the solution described in http://stackoverflow.com/a/88020
        # which is to use the parametric equation of an ellipse and to find
        # when dx/dt or dy/dt=0.

        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        tan_angle = np.tan(self.angle)

        t1 = np.arctan(-self.height * tan_angle / self.width)
        t2 = t1 + np.pi * u.rad

        dx1 = 0.5 * self.width * cos_angle * np.cos(t1) - 0.5 * self.height * sin_angle * np.sin(t1)
        dx2 = 0.5 * self.width * cos_angle * np.cos(t2) - 0.5 * self.height * sin_angle * np.sin(t2)

        if dx1 > dx2:
            dx1, dx2 = dx2, dx1

        t1 = np.arctan(self.height / tan_angle / self.width)
        t2 = t1 + np.pi * u.rad

        dy1 = 0.5 * self.height * cos_angle * np.sin(t1) + 0.5 * self.width * sin_angle * np.cos(t1)
        dy2 = 0.5 * self.height * cos_angle * np.sin(t2) + 0.5 * self.width * sin_angle * np.cos(t2)

        if dy1 > dy2:
            dy1, dy2 = dy2, dy1

        xmin = self.center.x + dx1
        xmax = self.center.x + dx2
        ymin = self.center.y + dy1
        ymax = self.center.y + dy2

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
            0.5 * self.width, 0.5 * self.height,
            self.angle.to(u.rad).value,
            use_exact, subpixels,
        )

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Ellipse`).
        """
        from matplotlib.patches import Ellipse
        xy = self.center.x, self.center.y
        width = self.width
        height = self.height
        # From the docstring: MPL expects "rotation in degrees (anti-clockwise)"
        angle = self.angle.to('deg').value
        return Ellipse(xy=xy, width=width, height=height, angle=angle, **kwargs)


class EllipseSkyRegion(SkyRegion):
    """
    An ellipse defined using sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the ellipse.
    width : `~astropy.units.Quantity`
        The width of the ellipse (before rotation) as an angle
    height : `~astropy.units.Quantity`
        The height of the ellipse (before rotation) as an angle
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the longitude axis
        of the celestial coordinates.
    """

    def __init__(self, center, width, height, angle=0. * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that height, width, and angle are angles
        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = [('width', self.width), ('height', self.height),
                             ('angle', self.angle)]

    def to_pixel(self, wcs):
        center, scale, north_angle = skycoord_to_pixel_scale_angle(self.center, wcs)
        # FIXME: The following line is needed to get a scalar PixCoord
        center = PixCoord(float(center.x), float(center.y))
        height = self.height.to('deg').value * scale
        width = self.width.to('deg').value * scale
        return EllipsePixelRegion(center, width, height,
                                  angle=self.angle + (north_angle - 90 * u.deg),
                                  meta=self.meta, visual=self.visual)
