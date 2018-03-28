# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import math
import operator

import numpy as np

from astropy import wcs
from astropy import coordinates
from astropy import units as u
from astropy.wcs.utils import pixel_to_skycoord

from .._utils.wcs_helpers import skycoord_to_pixel_scale_angle
from ..core import CompoundPixelRegion, CompoundSkyRegion, PixCoord
from ..shapes import CirclePixelRegion, CircleSkyRegion

__all__ = ['CircleAnnulusPixelRegion', 'CircleAnnulusSkyRegion']


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
        self._repr_params = [('inner radius', region1.radius),
                             ('outer radius', region2.radius)]

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
        return self.region2.bounding_box()

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
        self._repr_params = [('inner radius', region1.radius),
                             ('outer radius', region2.radius)]

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
