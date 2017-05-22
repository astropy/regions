# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy import units as u

from astropy.coordinates import Angle
from astropy.wcs.utils import pixel_to_skycoord

from ..core import PixCoord, PixelRegion, SkyRegion, Mask, BoundingBox
from .._geometry import rectangular_overlap_grid
from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle

__all__ = ['RectanglePixelRegion', 'RectangleSkyRegion']


class RectanglePixelRegion(PixelRegion):
    """
    A rectangle in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the rectangle.
    width : float
        The width of the rectangle (before rotation) in pixels
    height : float
        The height of the rectangle (before rotation) in pixels
    angle : `~astropy.units.Quantity`
        The rotation angle of the rectangle, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the x axis.

    Examples
    --------

    .. plot::
        :include-source:

        import numpy as np
        from astropy.coordinates import Angle
        from regions import PixCoord, RectanglePixelRegion
        import matplotlib.pyplot as plt

        x, y = 15, 10
        width, height = 8, 5
        angle = Angle(30, 'deg')

        fig, ax = plt.subplots(1, 1)

        center = PixCoord(x=x, y=y)
        reg = RectanglePixelRegion(center=center, width=width, height=height, angle=angle)
        patch = reg.as_patch(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(0, 30)
        plt.ylim(0, 20)
        ax.set_aspect('equal')
        plt.show()
    """

    def __init__(self, center, width, height, angle=0 * u.deg, meta=None, visual=None):
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
        return self.width * self.height

    def contains(self, pixcoord):
        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        dx = pixcoord.x - self.center.x
        dy = pixcoord.y - self.center.y
        dx_rot = cos_angle * dx + sin_angle * dy
        dy_rot = sin_angle * dx - cos_angle * dy
        return (np.abs(dx_rot) < self.width * 0.5) & (np.abs(dy_rot) < self.height * 0.5)

    def to_shapely(self):

        from shapely import affinity
        from shapely.geometry import Polygon

        x1 = self.center.x - self.width * 0.5
        y1 = self.center.y - self.height * 0.5
        x2 = self.center.x + self.width * 0.5
        y2 = self.center.y - self.height * 0.5
        x3 = self.center.x + self.width * 0.5
        y3 = self.center.y + self.height * 0.5
        x4 = self.center.x - self.width * 0.5
        y4 = self.center.y + self.height * 0.5

        rectangle = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)])

        return affinity.rotate(rectangle, self.angle.to(u.deg).value)

    def to_sky(self, wcs):
        # TODO: write a pixel_to_skycoord_scale_angle
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, north_angle = skycoord_to_pixel_scale_angle(center, wcs)
        width = Angle(self.width / scale, 'deg')
        height = Angle(self.height / scale, 'deg')
        return RectangleSkyRegion(center, width, height,
                                  angle=self.angle - (north_angle - 90 * u.deg),
                                  meta=self.meta, visual=self.visual)

    @property
    def bounding_box(self):
        """
        The minimal bounding box (`~regions.BoundingBox`) enclosing the
        exact rectangular region.
        """

        w2 = self.width / 2.
        h2 = self.height / 2.
        cos_angle = np.cos(self.angle)    # self.angle is a Quantity
        sin_angle = np.sin(self.angle)    # self.angle is a Quantity
        dx1 = abs(w2 * cos_angle - h2 * sin_angle)
        dy1 = abs(w2 * sin_angle + h2 * cos_angle)
        dx2 = abs(w2 * cos_angle + h2 * sin_angle)
        dy2 = abs(w2 * sin_angle - h2 * cos_angle)
        dx = max(dx1, dx2)
        dy = max(dy1, dy2)

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

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5 - self.center.x
        xmax = float(bbox.ixmax) - 0.5 - self.center.x
        ymin = float(bbox.iymin) - 0.5 - self.center.y
        ymax = float(bbox.iymax) - 0.5 - self.center.y

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = rectangular_overlap_grid(
            xmin, xmax, ymin, ymax, nx, ny,
            self.width, self.height,
            self.angle.to(u.rad).value,
            use_exact, subpixels,
        )

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Rectangle`).
        """
        from matplotlib.patches import Rectangle
        xy = self._lower_left_xy()
        width = self.width
        height = self.height
        # From the docstring: MPL expects "rotation in degrees (anti-clockwise)"
        angle = self.angle.to('deg').value
        return Rectangle(xy=xy, width=width, height=height, angle=angle, **kwargs)

    def _lower_left_xy(self):
        """
        Compute lower left `xy` position.

        This is used for the conversion to matplotlib in ``as_patch``

        Taken from http://photutils.readthedocs.io/en/latest/_modules/photutils/aperture/rectangle.html#RectangularAperture.plot
        """
        hw = self.width / 2.
        hh = self.height / 2.
        sint = np.sin(self.angle)
        cost = np.cos(self.angle)
        dx = (hh * sint) - (hw * cost)
        dy = -(hh * cost) - (hw * sint)
        x = self.center.x + dx
        y = self.center.y + dy
        return x, y


class RectangleSkyRegion(SkyRegion):
    """
    A rectangle in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangle.
    width : `~astropy.units.Quantity`
        The width of the rectangle (before rotation) as an angle
    height : `~astropy.units.Quantity`
        The height of the rectangle (before rotation) as an angle
    angle : `~astropy.units.Quantity`
        The rotation angle of the rectangle, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the longitude axis
        of the celestial coordinates.
    """

    def __init__(self, center, width, height, angle=0 * u.deg, meta=None, visual=None):
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
        width = self.width.to('deg').value * scale
        height = self.height.to('deg').value * scale
        return RectanglePixelRegion(center, width, height,
                                    angle=self.angle + (north_angle - 90 * u.deg),
                                    meta=self.meta, visual=self.visual)
