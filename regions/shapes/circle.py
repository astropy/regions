# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import math
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord
from ..core import PixelRegion, SkyRegion
from ..utils import skycoord_to_pixel_scale_angle
from ..utils.positions_to_extents import get_phot_extents
from ..utils.wcs_helpers import skycoord_to_pixel_scale_angle

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the circle.
    radius : float
        The radius of the circle
    """

    def __init__(self, center, radius, meta=None, visual=None):
        # TODO: test that center is a 0D PixCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
        return np.hypot(pixcoord.x - self.center.x,
                        pixcoord.y - self.center.y) < self.radius

    def to_shapely(self):
        return self.center.to_shapely().buffer(self.radius)

    def to_sky(self, wcs, mode='local', tolerance=None):
        if mode != 'local':
            raise NotImplementedError
        if tolerance is not None:
            raise NotImplementedError

        skypos = pixel_to_skycoord(self.center[0], self.center[1], wcs)
        xc, yc, scale, angle = skycoord_to_pixel_scale_angle(skypos, wcs)

        radius_sky = self.radius / scale
        return CircleSkyRegion(skypos, radius_sky)

    def to_mask(self, shape, mode='center'):
        # TODO: needs to be implemented

        positions = self.center
        radius = self.radius

        extents = np.zeros((len(positions), 4), dtype=int)

        extents[:, 0] = positions[:, 0] - radius + 0.5
        extents[:, 1] = positions[:, 0] + radius + 1.5
        extents[:, 2] = positions[:, 1] - radius + 0.5
        extents[:, 3] = positions[:, 1] + radius + 1.5

        ood_filter, extent, phot_extent = get_phot_extents(shape, positions,
                                                           extents)

        x_min, x_max, y_min, y_max = extent
        x_pmin, x_pmax, y_pmin, y_pmax = phot_extent

        if mode == 'center':
            use_exact = 0
            subpixels = 1
        elif mode == 'subpixel':
            use_exact = 0
        else:
            use_exact = 1
            subpixels = 1

        from .geometry import circular_overlap_grid

        for ind in range(len(positions)):
            fraction = circular_overlap_grid(x_pmin[ind], x_pmax[ind],
                                             y_pmin[ind], y_pmax[ind],
                                             x_max[ind] - x_min[ind],
                                             y_max[ind] - y_min[ind],
                                             radius, use_exact, subpixels)

        return fraction


    def as_patch(self, **kwargs):
        import matplotlib.patches as mpatches

        patch = mpatches.Circle(self.center, self.radius, **kwargs)
        return patch


class CircleSkyRegion(SkyRegion):
    """
    A circle in sky coordinates.

    Parameters
    ----------
    center : :class:`~astropy.coordinates.SkyCoord`
        The position of the center of the circle.
    radius : :class:`~astropy.units.Quantity`
        The radius of the circle in angular units
    """

    def __init__(self, center, radius, meta=None, visual=None):
        # TODO: test that center is a 0D SkyCoord
        self.center = center
        self.radius = radius
        self.meta = meta or {}
        self.visual = visual or {}

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, skycoord):
        return self.center.separation(skycoord) < self.radius

    def __repr__(self):
        clsnm = self.__class__.__name__
        coord = self.center
        rad = self.radius
        return '{clsnm}\nCenter:{coord}\nRadius:{rad}'.format(**locals())

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

        xc, yc, scale, angle = skycoord_to_pixel_scale_angle(self.center, wcs)
        radius_pix = (self.radius * scale)

        pixel_positions = np.array([xc, yc]).transpose()

        return CirclePixelRegion(pixel_positions, radius_pix)

    def as_patch(self, ax, **kwargs):
        import matplotlib.patches as mpatches

        val = self.center.icrs
        center = (val.ra.to('deg').value, val.dec.to('deg').value)

        temp = dict(transform=ax.get_transform('icrs'),
                    radius=self.radius.to('deg').value)
        kwargs.update(temp)
        for key, value in self.visual.items():
            kwargs.setdefault(key, value)
        patch = mpatches.Circle(center, **kwargs)

        return patch
