# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy import units as u

from ..core import PixelRegion, SkyRegion, Mask
from .._geometry import rectangular_overlap_grid

__all__ = ['RectanglePixelRegion', 'RectangleSkyRegion']


class RectanglePixelRegion(PixelRegion):
    """
    A rectangle in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the rectangle.
    height : float
        The height of the rectangle
    width : float
        The width of the rectangle
    angle : `~astropy.units.Quantity`
        The rotation of the rectangle. If set to zero (the default), the width
        is lined up with the x axis.
    """

    def __init__(self, center, height, width, angle=0 * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that angle is an angle
        self.center = center
        self.height = height
        self.width = width
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            height=self.height,
            width=self.width,
            angle=self.angle,
        )
        fmt = '{name}\ncenter: {center}\nheight: {height}\nwidth: {width}\nangle: {angle}'
        return fmt.format(**data)

    @property
    def area(self):
        return self.width * self.height

    def contains(self, pixcoord):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_shapely(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_sky(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_mask(self, mode='center', subpixels=5):

        # NOTE: assumes this class represents a single circle

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find exact bounds
        # FIXME: this is not the minimal bounding box, and can be optimized
        radius = np.hypot(self.width / 2, self.height / 2)
        xmin = self.center.x - radius
        xmax = self.center.x + radius
        ymin = self.center.y - radius
        ymax = self.center.y + radius

        # Find range of pixels
        imin = int(xmin)
        imax = int(xmax) + 1
        jmin = int(ymin)
        jmax = int(ymax) + 1

        # Find mask size
        nx = imax - imin + 1
        ny = jmax - jmin + 1

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(imin) - 0.5 - self.center.x
        xmax = float(imax) + 0.5 - self.center.x
        ymin = float(jmin) - 0.5 - self.center.y
        ymax = float(jmax) + 0.5 - self.center.y

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = rectangular_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                            self.width, self.height,
                                            self.angle.to(u.deg).value,
                                            use_exact, subpixels)

        return Mask(fraction, origin=(jmin, imin))

    def as_patch(self, **kwargs):
        # TODO: needs to be implemented
        raise NotImplementedError


class RectangleSkyRegion(SkyRegion):
    """
    A rectangle in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangle.
    height : `~astropy.units.Quantity`
        The height radius of the rectangle
    width : `~astropy.units.Quantity`
        The width radius of the rectangle
    angle : `~astropy.units.Quantity`
        The rotation of the rectangle. If set to zero (the default), the width
        is lined up with the longitude axis of the celestial coordinates.
    """

    def __init__(self, center, height, width, angle=0 * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that height, width, and angle are angles
        self.center = center
        self.height = height
        self.width = width
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            height=self.height,
            width=self.width,
            angle=self.angle,
        )
        fmt = '{name}\ncenter: {center}\nheight: {height}\nwidth: {width}\nangle: {angle}'
        return fmt.format(**data)

    @property
    def area(self):
        return self.width * self.height

    def contains(self, skycoord):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_pixel(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_patch(self, **kwargs):
        # TODO: needs to be implemented
        raise NotImplementedError
