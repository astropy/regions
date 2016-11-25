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
    minor : float
        Minor radius
    major : float
        Major radius
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse.
        If set to zero (the default), the major
        axis is lined up with the x axis.
    """

    def __init__(self, center, minor, major, angle=0. * u.deg, meta=None,
                 visual=None):
        # TODO: use quantity_input to check that angle is an angle
        self.center = center
        self.minor = minor
        self.major = major
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            minor=self.minor,
            major=self.major,
            angle=self.angle,
        )
        fmt = '{name}\ncenter: {center}\nminor: {minor}\nmajor: {major}\nangle: {angle}'
        return fmt.format(**data)

    @property
    def area(self):
        return math.pi * self.minor * self.major

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

        # Find exact bounds
        # FIXME: this is not the minimal bounding box, and can be optimized
        xmin = self.center.x - max(self.major, self.minor)
        xmax = self.center.x + max(self.major, self.minor)
        ymin = self.center.y - max(self.major, self.minor)
        ymax = self.center.y + max(self.major, self.minor)

        # Find range of pixels. We use round here because if the region extends
        # to e.g. -0.4, it's enough to set the bounding box lower value to 0
        # because the 0-th pixel goes from -0.5 to 0.5. At the upper end we add
        # 1 because the upper limits need to be exlcusive.
        ixmin = round(xmin)
        ixmax = round(xmax) + 1
        iymin = round(ymin)
        iymax = round(ymax) + 1

        return BoundingBox(ixmin, ixmax, iymin, iymax)

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

        fraction = elliptical_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                           self.major, self.minor,
                                           self.angle.to(u.deg).value,
                                           use_exact, subpixels)

        return Mask(fraction, bbox=bbox)

    def as_patch(self, **kwargs):
        # TODO: needs to be implemented
        raise NotImplementedError


class EllipseSkyRegion(SkyRegion):
    """
    An ellipse in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Center position
    minor : `~astropy.units.Quantity`
        Minor radius
    major : `~astropy.units.Quantity`
        Major radius
    angle : `~astropy.units.Quantity`
        The rotation angle of the ellipse.
        If set to zero (the default), the major
        axis is lined up with the longitude axis of the celestial coordinates.
    """

    def __init__(self, center, minor, major, angle=0. * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that height, width, and angle are angles
        self.center = center
        self.minor = minor
        self.major = major
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            minor=self.minor,
            major=self.major,
            angle=self.angle,
        )
        fmt = '{name}\ncenter: {center}\nminor: {minor}\nmajor: {major}\nangle: {angle}'
        return fmt.format(**data)

    @property
    def area(self):
        # TODO: needs to be implemented
        raise NotImplementedError

    def contains(self, skycoord):
        # TODO: needs to be implemented
        raise NotImplementedError

    def to_pixel(self, wcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_patch(self, **kwargs):
        # TODO: needs to be implemented
        raise NotImplementedError
