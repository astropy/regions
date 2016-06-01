# Licensed under a 3-clause BSD style license - see LICENSE.rst

import math

import numpy as np

from astropy import wcs
from astropy import coordinates
from astropy import units as u

from ..core import PixelRegion, SkyRegion
from ..utils.wcs_helpers import skycoord_to_pixel_scale_angle


class AnnulusPixelRegion(PixelRegion):
    """
    An annulus in pixel coordinates.

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
        # TODO: test that center is a 0D PixCoord
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or {}
        self.visual = visual or {}

    @property
    def area(self):
        return math.pi * (self.outer_radius ** 2 - self.inner_radius ** 2)

    def contains(self, pixcoord):
        distance = np.hypot(pixcoord.x - self.center.x, pixcoord.y - self.center.y) 
        return distance < self.outer_radius and distance > self.inner_radius

    def to_shapely(self):
        # TODO: needs to be implemented
        raise NotImplementedError("")

    def to_sky(self, mywcs, mode='local', tolerance=None):
        # TODO: needs to be implemented
        raise NotImplementedError("")

    def to_mask(self, mode='center'):
        # TODO: needs to be implemented
        raise NotImplementedError("")

    def as_patch(self, **kwargs):
        import matplotlib.patches as mpatches
        # TODO: needs to be implemented
        raise NotImplementedError("")


class AnnulusSkyRegion(SkyRegion):
    """
    An annulus in sky coordinates.

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
        # TODO: test that center is a 0D SkyCoord
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or {}
        self.visual = visual or {}

    @property
    def area(self):
        return math.pi * (self.outer_radius ** 2 - self.inner_radius ** 2)

    def contains(self, skycoord):
        distance = self.center.separation(skycoord)
        return  distance < self.outer_radius and distance > self.inner_radius

    def __repr__(self):
        clsnm = self.__class__.__name__
        coord = self.center
        irad = self.inner_radius
        orad = self.outer_radius
        ss = '{clsnm}\nCenter:{coord}\nInner Radius:{irad}\nOuter Radius:{orad}'
        return ss.format(**locals())

    def to_pixel(self, mywcs, mode='local', tolerance=None):
        """
        Given a WCS, convert the circle to a best-approximation circle in pixel
        dimensions.

        Parameters
        ----------
        mywcs : `~astropy.wcs.WCS`
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
            raise NotImplementedError()
        if tolerance is not None:
            raise NotImplementedError()

        # TODO: needs to be implemented
        raise NotImplementedError("")

    def as_patch(self, ax, **kwargs):
        import matplotlib.patches as mpatches
        # TODO: needs to be implemented
        raise NotImplementedError("")
