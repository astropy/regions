# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from ..core import PixCoord, PixelRegion, SkyRegion, BoundingBox
from ..core.attributes import ScalarPix, ScalarSky, RegionMeta, RegionVisual

__all__ = ['PointPixelRegion', 'PointSkyRegion']


class PointPixelRegion(PixelRegion):
    """
    A point position in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the point
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

    Examples
    --------

    .. plot::
        :include-source:

        from regions import PixCoord, PointPixelRegion, RegionVisual
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)
        regs = []
        regs.append(PointPixelRegion(PixCoord(2, 2), visual=RegionVisual(symbol='D')))
        regs.append(PointPixelRegion(PixCoord(2, 3), visual=RegionVisual(symbol='*')))
        regs.append(PointPixelRegion(PixCoord(3, 3), visual=RegionVisual(symbol='^')))
        regs.append(PointPixelRegion(PixCoord(3, 2), visual=RegionVisual(symbol='*')))
        regs.append(PointPixelRegion(PixCoord(2, 4), visual=RegionVisual(symbol='x')))
        regs.append(PointPixelRegion(PixCoord(4, 2)))

        for reg in regs:
            reg.plot(ax=ax)

        plt.xlim(0, 6)
        plt.ylim(0, 6)
        ax.set_aspect('equal')

    """

    center = ScalarPix('center')

    def __init__(self, center, meta=None, visual=None):
        self.center = center
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = tuple()

    @property
    def area(self):
        return 0

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            in_reg = False
        else:
            in_reg = np.zeros(pixcoord.x.shape, dtype=bool)

        if self.meta.get('include', True):
            # in_reg = False, always.  Points do not include anything
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return PointSkyRegion(center)

    @property
    def bounding_box(self):
        return BoundingBox.from_float(self.center.x, self.center.x,
                                      self.center.y, self.center.y)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Matplotlib Line2D object for this region (`matplotlib.lines.Line2D`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.lines.Line2D` object accepts

        Returns
        -------
        point : `~matplotlib.lines.Line2D`
            Matplotlib Line2D object.
        """

        from matplotlib.lines import Line2D

        mpl_params = self.mpl_properties_default('LINE2D')
        mpl_params.update(kwargs)

        point = Line2D([self.center.x - origin[0]], [self.center.y - origin[1]],
                       **mpl_params)

        return point

    def rotate(self, center, angle):
        """Make a rotated region.

        Rotates counter-clockwise for positive ``angle``.

        Parameters
        ----------
        center : `PixCoord`
            Rotation center point
        angle : `~astropy.coordinates.Angle`
            Rotation angle

        Returns
        -------
        region : `PointPixelRegion`
            Rotated region (an independent copy)
        """
        center = self.center.rotate(center, angle)
        return self.copy(center=center)


class PointSkyRegion(SkyRegion):
    """
    A pixel region in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the point
    meta : `regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """

    center = ScalarSky('center')

    def __init__(self, center, meta=None, visual=None):
        self.center = center
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        self._repr_params = tuple()

    def contains(self, skycoord, wcs):
        if self.meta.get('include', True):
            # points never include anything
            return False
        else:
            return True

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return PointPixelRegion(center)
