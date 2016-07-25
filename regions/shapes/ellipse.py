# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import math
from astropy import units as u
from ..core import PixelRegion, SkyRegion

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

    def to_mask(self, mode='center'):
        # TODO: needs to be implemented
        raise NotImplementedError

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
