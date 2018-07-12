# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import operator
import six
import abc


from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord

from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle
from ..core import CompoundPixelRegion, CompoundSkyRegion, PixCoord
from ..shapes.circle import CirclePixelRegion, CircleSkyRegion
from ..shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from ..shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion

__all__ = ['CircleAnnulusPixelRegion', 'CircleAnnulusSkyRegion',
           'EllipseAnnulusPixelRegion', 'EllipseAnnulusSkyRegion',
           'RectangleAnnulusPixelRegion', 'RectangleAnnulusSkyRegion'
           ]


class CircleAnnulusPixelRegion(CompoundPixelRegion):
    """
    A circular annulus in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the annulus.
    inner_radius : float
        The inner radius of the annulus
    outer_radius : float
        The outer radius of the annulus
    """

    def __init__(self, center, inner_radius, outer_radius, meta=None, visual=None):
        if inner_radius > outer_radius:
            raise ValueError('Outer radius should be larger than inner radius.')
        region1 = CirclePixelRegion(center, inner_radius)
        region2 = CirclePixelRegion(center, outer_radius)
        super(CircleAnnulusPixelRegion, self).__init__(
            region1=region1, region2=region2, operator=operator.xor, meta=meta, visual=visual)
        self._repr_params = ('inner_radius', 'outer_radius')

    @property
    def center(self):
        return self.region1.center

    @property
    def inner_radius(self):
        return self.region1.radius

    @property
    def outer_radius(self):
        return self.region2.radius

    @property
    def area(self):
        return self.region2.area - self.region1.area

    @property
    def bounding_box(self):
        return self.region2.bounding_box

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, _ = skycoord_to_pixel_scale_angle(center, wcs)
        inner_radius = self.inner_radius / scale * u.deg
        outer_radius = self.outer_radius / scale * u.deg
        return CircleAnnulusSkyRegion(center, inner_radius, outer_radius, self.meta, self.visual)


class CircleAnnulusSkyRegion(CompoundSkyRegion):
    """
    A circular annulus in sky coordinates.

    Parameters
    ----------
    center : :class:`~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_radius : :class:`~astropy.units.Quantity`
        The inner radius of the annulus in angular units
    outer_radius : :class:`~astropy.units.Quantity`
        The outer radius of the annulus in angular units
    """

    def __init__(self, center, inner_radius, outer_radius, meta=None, visual=None):
        if inner_radius > outer_radius:
            raise ValueError('Outer radius should be larger than inner radius.')
        region1 = CircleSkyRegion(center, inner_radius)
        region2 = CircleSkyRegion(center, outer_radius)
        super(CircleAnnulusSkyRegion, self).__init__(
            region1=region1, operator=operator.xor, region2=region2, meta=meta, visual=visual)
        self._repr_params = ('inner_radius', 'outer_radius')

    @property
    def center(self):
        return self.region1.center

    @property
    def inner_radius(self):
        return self.region1.radius

    @property
    def outer_radius(self):
        return self.region2.radius

    def to_pixel(self, wcs):
        center, scale, _ = skycoord_to_pixel_scale_angle(self.center, wcs)
        # FIXME: The following line is needed to get a scalar PixCoord
        center = PixCoord(float(center.x), float(center.y))
        inner_radius = self.inner_radius.to('deg').value * scale
        outer_radius = self.outer_radius.to('deg').value * scale
        return CircleAnnulusPixelRegion(center, inner_radius, outer_radius, self.meta, self.visual)


@six.add_metaclass(abc.ABCMeta)
class AsymmetricAnnulusPixelRegion(CompoundPixelRegion):
    """
    An asymmetrical annulus in pixel coordinates.

    region1: `~regions.PixelRegion` object
        The inner asymmetric region
    region2: `~regions.PixelRegion` object
        The outer asymmetric region
    meta: `region.RegionMeta` object
        The meta attributes of the region.
    visual: `region.RegionVisual` object
        The visual meta attributes of the region.
    """

    def __init__(self, region1, region2, meta=None, visual=None):

        super(AsymmetricAnnulusPixelRegion, self).__init__(
            region1=region1, region2=region2, operator=operator.xor, meta=meta, visual=visual)

        self._repr_params = ('inner_width', 'inner_height',
                             'outer_width', 'outer_height', 'angle')

    @property
    def center(self):
        return self.region1.center

    @property
    def inner_width(self):
        return self.region1.width

    @property
    def inner_height(self):
        return self.region1.height

    @property
    def outer_width(self):
        return self.region2.width

    @property
    def outer_height(self):
        return self.region2.height

    @property
    def angle(self):
        return self.region2.angle

    @property
    def area(self):
        return self.region2.area - self.region1.area

    @property
    def bounding_box(self):
        return self.region2.bounding_box

    def to_sky_args(self, wcs):

        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, scale, north_angle = skycoord_to_pixel_scale_angle(center, wcs)

        inner_width = self.inner_width / scale * u.deg
        inner_height = self.inner_height / scale * u.deg
        outer_width = self.outer_width / scale * u.deg
        outer_height = self.outer_height / scale * u.deg
        angle = self.angle - (north_angle - 90 * u.deg)

        return center, inner_width, inner_height, outer_width, outer_height, angle


@six.add_metaclass(abc.ABCMeta)
class AsymmetricAnnulusSkyRegion(CompoundSkyRegion):
    """
    A rectangular annulus in `~astropy.coordinates.SkyCoord` coordinates.

    Parameters
    ----------
    region1: `~regions.SkyRegion` object
        The inner asymmetric region
    region2: `~regions.SkyRegion` object
        The outer asymmetric region
    meta: `region.RegionMeta` object
        The meta attributes of the region.
    visual: `region.RegionVisual` object
        The visual meta attributes of the region.
    """

    def __init__(self, region1, region2, meta=None, visual=None):

        super(AsymmetricAnnulusSkyRegion, self).__init__(region1=region1, region2=region2, operator=operator.xor,
                meta=meta, visual=visual)

        self._repr_params = ('inner_width', 'inner_height',
                             'outer_width', 'outer_height', 'angle')

    @property
    def center(self):
        return self.region1.center

    @property
    def inner_width(self):
        return self.region1.width

    @property
    def inner_height(self):
        return self.region1.height

    @property
    def outer_width(self):
        return self.region2.width

    @property
    def outer_height(self):
        return self.region2.height

    @property
    def angle(self):
        return self.region2.angle

    def to_pixel_args(self, wcs):
        center, scale, north_angle = skycoord_to_pixel_scale_angle(self.center, wcs)
        # FIXME: The following line is needed to get a scalar PixCoord
        center = PixCoord(float(center.x), float(center.y))

        inner_width = self.inner_width.to('deg').value * scale
        inner_height = self.inner_height.to('deg').value * scale
        outer_width = self.outer_width.to('deg').value * scale
        outer_height = self.outer_height.to('deg').value * scale
        angle = self.angle + (north_angle - 90 * u.deg)

        return center,inner_width, inner_height, outer_width, outer_height, angle


class EllipseAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A elliptical annulus in pixel coordinates.

    Parameters
    ----------
    center: `~regions.PixCoord`
        The position of the center of the elliptical annulus.
    inner_width: float
        The inner width of the elliptical annulus (before rotation) in pixels
    inner_height: float
        The inner height of the elliptical annulus (before rotation) in pixels
    outer_width: float
        The outer width of the elliptical annulus (before rotation) in pixels
    outer_height: float
        The outer height of the elliptical annulus (before rotation) in pixels
    angle: `~astropy.units.Quantity`
        The rotation angle of the elliptical annulus, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the x axis.
    """

    def __init__(self, center, inner_width, inner_height, outer_width,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        if inner_width > outer_width:
            raise ValueError('Outer width should be larger than inner width.')

        if inner_height > outer_height:
            raise ValueError('Outer height should be larger than inner height.')

        region1 = EllipsePixelRegion(center, inner_width, inner_height, angle)
        region2 = EllipsePixelRegion(center, outer_width, outer_height, angle)

        super(EllipseAnnulusPixelRegion, self).__init__(
            region1=region1, region2=region2, meta=meta, visual=visual)

    def to_sky(self, wcs):
        return EllipseAnnulusSkyRegion(*self.to_sky_args(wcs),
                                       meta=self.meta, visual=self.visual)


class EllipseAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A elliptical annulus in `~astropy.coordinates.SkyCoord` coordinates.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The position of the center of the elliptical annulus.
    inner_width: `~astropy.units.Quantity`
        The inner width of the elliptical annulus (before rotation) as angle
    inner_height: `~astropy.units.Quantity`
        The inner height of the elliptical annulus (before rotation) as angle
    outer_width: `~astropy.units.Quantity`
        The outer width of the elliptical annulus (before rotation) as angle
    outer_height: `~astropy.units.Quantity`
        The outer height of the elliptical annulus (before rotation) as angle
    angle: `~astropy.units.Quantity`, optional
        The rotation angle of the elliptical annulus, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the longitude axis
        of the celestial coordinates
    """

    def __init__(self, center, inner_width, inner_height, outer_width,
                    outer_height, angle=0 * u.deg, meta=None, visual=None):
        if inner_width > outer_width:
            raise ValueError('Outer width should be larger than inner width.')

        if inner_height > outer_height:
            raise ValueError('Outer height should be larger than inner height.')

        region1 = EllipseSkyRegion(center, inner_width, inner_height, angle)
        region2 = EllipseSkyRegion(center, outer_width, outer_height, angle)

        super(EllipseAnnulusSkyRegion, self).__init__(
            region1=region1, region2=region2, meta=meta, visual=visual)

    def to_pixel(self, wcs):
        return EllipseAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                         meta=self.meta, visual=self.visual)


class RectangleAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A rectangular annulus in pixel coordinates.

    Parameters
    ----------
    center: `~regions.PixCoord`
        The position of the center of the rectangular annulus.
    inner_width: float
        The inner width of the rectangular annulus (before rotation) in pixels
    inner_height: float
        The inner height of the rectangular annulus (before rotation) in pixels
    outer_width: float
        The outer width of the rectangular annulus (before rotation) in pixels
    outer_height: float
        The outer height of the rectangular annulus (before rotation) in pixels
    angle: `~astropy.units.Quantity`
        The rotation angle of the rectangular annulus, measured anti-clockwise.
        If set to zero (the default), the width axis is lined up with the x axis.
    """

    def __init__(self, center, inner_width, inner_height, outer_width,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        if inner_width > outer_width:
            raise ValueError('Outer width should be larger than inner width.')

        if inner_height > outer_height:
            raise ValueError('Outer height should be larger than inner height.')

        region1 = RectanglePixelRegion(center, inner_width, inner_height, angle)
        region2 = RectanglePixelRegion(center, outer_width, outer_height, angle)

        super(RectangleAnnulusPixelRegion, self).__init__(
            region1=region1, region2=region2, meta=meta, visual=visual)

    def to_sky(self, wcs):
        return RectangleAnnulusSkyRegion(*self.to_sky_args(wcs),
                                         meta=self.meta, visual=self.visual)


class RectangleAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A rectangular annulus in `~astropy.coordinates.SkyCoord` coordinates.

    Parameters
    ----------
    center: `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangular annulus.
    inner_width: `~astropy.units.Quantity`
        The inner width of the rectangular annulus (before rotation) as angle
    inner_height: `~astropy.units.Quantity`
        The inner height of the rectangular annulus (before rotation) as angle
    outer_width: `~astropy.units.Quantity`
        The outer width of the rectangular annulus (before rotation) as angle
    outer_height: `~astropy.units.Quantity`
        The outer height of the rectangular annulus (before rotation) as angle
    angle: `~astropy.units.Quantity`, optional
        The rotation angle of the rectangular annulus, measured anti-clockwise. If set to
        zero (the default), the width axis is lined up with the longitude axis
        of the celestial coordinates
    """

    def __init__(self, center, inner_width, inner_height, outer_width,
                    outer_height, angle=0 * u.deg, meta=None, visual=None):
        if inner_width > outer_width:
            raise ValueError('Outer width should be larger than inner width.')

        if inner_height > outer_height:
            raise ValueError('Outer height should be larger than inner height.')

        region1 = RectangleSkyRegion(center, inner_width, inner_height, angle)
        region2 = RectangleSkyRegion(center, outer_width, outer_height, angle)

        super(RectangleAnnulusSkyRegion, self).__init__(
            region1=region1, region2=region2, meta=meta, visual=visual)

    def to_pixel(self, wcs):
        return RectangleAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                           meta=self.meta, visual=self.visual)
