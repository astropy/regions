# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import math
from astropy.coordinates import Angle
from astropy.wcs.utils import pixel_to_skycoord
from astropy import units as u
from ..core import PixelRegion, SkyRegion
from ..utils.wcs_helpers import skycoord_to_pixel_scale_angle

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        Center position
    radius : float
        Radius
    """

    def __init__(self, center, radius, meta=None, visual=None):
        # TODO: test that center is a 0D PixCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            radius=self.radius,
        )
        fmt = '{name}\ncenter: {center}\nradius: {radius}'
        return fmt.format(**data)

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
        return self.center.separation(pixcoord) < self.radius

    def to_shapely(self):
        return self.center.to_shapely().buffer(self.radius)

    def to_sky(self, wcs, mode='local', tolerance=None):
        if mode != 'local':
            raise NotImplementedError
        if tolerance is not None:
            raise NotImplementedError

        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        # TODO: this is just called to compute `scale`
        # This is inefficient ... we should have that as a separate function.
        _, scale, _ = skycoord_to_pixel_scale_angle(center, wcs)

        radius = Angle(self.radius / scale, 'deg')
        return CircleSkyRegion(center, radius)

    def to_mask(self, mode='center'):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_patch(self, **kwargs):

        import matplotlib.patches as mpatches

        xy = self.center.x, self.center.y
        radius = self.radius
        patch = mpatches.Circle(xy=xy, radius=radius, **kwargs)
        return patch


class CircleSkyRegion(SkyRegion):
    """
    A circle in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        Center position
    radius : `~astropy.units.Quantity`
        Radius in angular units
    """

    def __init__(self, center, radius, meta=None, visual=None):
        # TODO: test that center is a 0D SkyCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}

    def __repr__(self):
        data = dict(
            name=self.__class__.__name__,
            center=self.center,
            radius=self.radius,
        )
        fmt = '{name}\ncenter: {center}\nradius: {radius}'
        return fmt.format(**data)

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, skycoord):
        return self.center.separation(skycoord) < self.radius

    def to_pixel(self, wcs, mode='local', tolerance=None):
        """
        Given a WCS, convert the circle to a best-approximation circle in pixel
        dimensions.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
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
            raise NotImplementedError
        if tolerance is not None:
            raise NotImplementedError

        center, scale, _ = skycoord_to_pixel_scale_angle(self.center, wcs)
        radius = self.radius.to('deg').value * scale

        return CirclePixelRegion(center, radius)

    def as_patch(self, ax, **kwargs):

        from wcsaxes.patches import SphericalCircle

        c_icrs = self.center.icrs

        for key, value in self.visual.items():
            kwargs.setdefault(key, value)

        patch = SphericalCircle((c_icrs.ra, c_icrs.dec), self.radius,
                                vertex_unit=u.degree,
                                transform=ax.get_transform('icrs'), **kwargs)

        return patch
