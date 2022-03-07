# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines point regions in both pixel and sky coordinates.
"""

from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import numpy as np

from ..core.attributes import ScalarPixCoord, ScalarSkyCoord
from ..core.bounding_box import RegionBoundingBox
from ..core.core import PixelRegion, SkyRegion
from ..core.metadata import RegionMeta, RegionVisual
from ..core.pixcoord import PixCoord


__all__ = ['PointPixelRegion', 'PointSkyRegion']


class PointPixelRegion(PixelRegion):
    """
    A point position in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the point.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, PointPixelRegion, RegionVisual
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        regs = []
        regs.append(PointPixelRegion(PixCoord(2, 2),
                    visual=RegionVisual(marker='D')))
        regs.append(PointPixelRegion(PixCoord(2, 3),
                    visual=RegionVisual(marker='+')))
        regs.append(PointPixelRegion(PixCoord(3, 3),
                    visual=RegionVisual(marker='^')))
        regs.append(PointPixelRegion(PixCoord(3, 2),
                    visual=RegionVisual(marker='*')))
        regs.append(PointPixelRegion(PixCoord(2, 4),
                    visual=RegionVisual(marker='x')))
        regs.append(PointPixelRegion(PixCoord(4, 2)))
        for reg in regs:
            reg.plot(ax=ax)

        ax.set_xlim(0, 6)
        ax.set_ylim(0, 6)
        ax.set_aspect('equal')
    """

    _params = ('center',)
    _mpl_artist = 'Line2D'
    center = ScalarPixCoord(description='The point pixel position.')

    def __init__(self, center, meta=None, visual=None):
        self.center = center
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
        return 0.0

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            in_reg = False
        else:
            in_reg = np.zeros(pixcoord.x.shape, dtype=bool)

        if self.meta.get('include', True):
            # in_reg = False, always.  Points do not include anything.
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return PointSkyRegion(center, meta=self.meta.copy(),
                              visual=self.visual.copy())

    @property
    def bounding_box(self):
        return RegionBoundingBox.from_float(self.center.x, self.center.x,
                                            self.center.y, self.center.y)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib Line2D object for this region
        (`matplotlib.lines.Line2D`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : `dict`
            Any keyword arguments accepted by
            `~matplotlib.lines.Line2D`. These keywords will override any
            visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.lines.Line2D`
            A matplotlib Line2D object.
        """
        from matplotlib.lines import Line2D

        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Line2D([self.center.x - origin[0]],
                      [self.center.y - origin[1]], **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `PointPixelRegion`
            The rotated region (which is an independent copy).
        """
        center = self.center.rotate(center, angle)
        return self.copy(center=center)


class PointSkyRegion(SkyRegion):
    """
    A pixel region in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the point.
    meta : `regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('center',)
    center = ScalarSkyCoord(description=('The point position as a sky '
                                         'coordinate.'))

    def __init__(self, center, meta=None, visual=None):
        self.center = center
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def contains(self, skycoord, wcs):  # pylint: disable=unused-argument
        if self.meta.get('include', True):
            # points never include anything
            return False
        else:
            return True

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return PointPixelRegion(center, meta=self.meta.copy(),
                                visual=self.visual.copy())
