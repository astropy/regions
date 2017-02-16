# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import math
from astropy import units as u
from ..core import PixelRegion, SkyRegion, Mask, BoundingBox
from .._geometry import elliptical_overlap_grid

__all__ = ['EllipsePixelRegion', 'EllipseSkyRegion']


class EllipsePixelRegion(PixelRegion):
    """
    An ellipse in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        Center position
    major : float
        Minor radius
    minor : float
        Major radius
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
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_shapely(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_sky(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError

    @property
    def bounding_box(self):
        """
        Bounding box (`~regions.BoundingBox`).
        """
        # Find exact bounds
        # FIXME: this is not the minimal bounding box, and can be optimized
        xmin = self.center.x - max(self.major, self.minor)
        xmax = self.center.x + max(self.major, self.minor)
        ymin = self.center.y - max(self.major, self.minor)
        ymax = self.center.y + max(self.major, self.minor)

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
            self.angle.to(u.deg).value,
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
    An ellipse in sky coordinates.

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

    @property
    def area(self):
        """Region sky area approximation (`~astropy.units.Quantity`)"""
        return math.pi * self.major * self.minor

    def contains(self, skycoord):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_pixel(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError
