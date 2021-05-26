# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines line regions in both pixel and sky coordinates.
"""

from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import numpy as np

from ..core.attributes import ScalarPix, ScalarSky
from ..core.bounding_box import BoundingBox
from ..core.core import PixelRegion, SkyRegion
from ..core.metadata import RegionMeta, RegionVisual
from ..core.pixcoord import PixCoord

__all__ = ['LinePixelRegion', 'LineSkyRegion']


class LinePixelRegion(PixelRegion):
    """
    A line in pixel coordinates.

    Parameters
    ----------
    start : `~regions.PixCoord`
        The start position.
    end : `~regions.PixCoord`
        The end position.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, LinePixelRegion
        import matplotlib.pyplot as plt

        x1, y1 = 15, 10
        x2, y2 = 20, 25

        fig, ax = plt.subplots(1, 1)

        start = PixCoord(x=x1, y=y1)
        end = PixCoord(x=x2, y=y2)
        reg = LinePixelRegion(start=start, end=end)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(0, 30)
        plt.ylim(0, 30)
        ax.set_aspect('equal')
    """

    _params = ('start', 'end')
    start = ScalarPix('start')
    end = ScalarPix('end')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
        return 0

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            in_reg = False
        else:
            in_reg = np.zeros(pixcoord.x.shape, dtype=bool)

        if self.meta.get('include', True):
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_sky(self, wcs):
        start = pixel_to_skycoord(self.start.x, self.start.y, wcs)
        end = pixel_to_skycoord(self.end.x, self.end.y, wcs)
        return LineSkyRegion(start, end)

    @property
    def bounding_box(self):
        xmin = min(self.start.x, self.end.x)
        xmax = max(self.start.x, self.end.x)
        ymin = min(self.start.y, self.end.y)
        ymax = max(self.start.y, self.end.y)

        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Arrow`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Arrow`. These keywords will override
            any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Arrow`
            A matplotlib line patch.
        """
        # Note: Long term we want to support DS9 lines with arrow heads.
        # We may want to use Line2D instead of arrow for lines because the
        # width of the arrow is non-scalable in patches.
        from matplotlib.patches import Arrow

        x = self.start.x - origin[0]
        y = self.start.y - origin[1]
        dx = self.end.x - self.start.x
        dy = self.end.y - self.start.y
        kwargs.setdefault("width", 0.1)

        mpl_kwargs = self._define_mpl_kwargs(artist='Patch')
        mpl_kwargs.update(kwargs)

        return Arrow(x, y, dx, dy, **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Postive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `LinePixelRegion`
            The rotated region (which is an independent copy).
        """
        start = self.start.rotate(center, angle)
        end = self.end.rotate(center, angle)
        return self.copy(start=start, end=end)


class LineSkyRegion(SkyRegion):
    """
    A line in sky coordinates.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        The start position.
    end : `~astropy.coordinates.SkyCoord`
        The end position.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('start', 'end')
    start = ScalarSky('start')
    end = ScalarSky('end')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def contains(self, skycoord, wcs):  # pylint: disable=unused-argument
        if self.meta.get('include', True):
            # lines never contain anything
            return False
        else:
            return True

    def to_pixel(self, wcs):
        start_x, start_y = skycoord_to_pixel(self.start, wcs=wcs)
        start = PixCoord(start_x, start_y)
        end_x, end_y = skycoord_to_pixel(self.end, wcs=wcs)
        end = PixCoord(end_x, end_y)
        return LinePixelRegion(start, end)
