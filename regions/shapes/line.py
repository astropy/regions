# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import astropy.units as u
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from ..core import PixCoord, PixelRegion, SkyRegion, BoundingBox
from ..core.attributes import ScalarPix, ScalarSky


__all__ = ['LinePixelRegion', 'LineSkyRegion']


class LinePixelRegion(PixelRegion):
    """
    A line in pixel coordinates.

    Parameters
    ----------
    start : `~regions.PixCoord`
        Start position
    end : `~regions.PixCoord`
        End position
    """

    start = ScalarPix('start')
    end = ScalarPix('end')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = ('start', 'end')

    @property
    def area(self):
        """Region area (float)."""
        return 0 * u.sr

    def contains(self, pixcoord):
        if pixcoord.isscalar:
            return False
        else:
            return np.zeros(pixcoord.x.shape, dtype=bool)

    def to_shapely(self):
        from shapely.geometry import LineString
        return LineString([(self.start.x, self.start.y), (self.end.x, self.end.y)])

    def to_sky(self, wcs, mode='local', tolerance=None):
        start = pixel_to_skycoord(self.start.x, self.start.y, wcs)
        end = pixel_to_skycoord(self.end.x, self.end.y, wcs)
        return LineSkyRegion(start, end)

    @property
    def bounding_box(self):
        """Bounding box (`~regions.BoundingBox`)."""
        xmin = min(self.start.x, self.end.x)
        xmax = max(self.start.x, self.end.x)
        ymin = min(self.start.y, self.end.y)
        ymax = max(self.start.y, self.end.y)

        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_patch(self, **kwargs):
        """Matplotlib patch object for this region (`matplotlib.patches.Line`)."""
        # Long term we want to support DS9 lines with arrow heads
        from matplotlib.patches import Arrow
        x = self.start.x
        y = self.start.y
        dx = self.end.x - self.start.x
        dy = self.end.y - self.start.y
        return Arrow(x, y, dx, dy, **kwargs)


class LineSkyRegion(SkyRegion):
    """
    A line in sky coordinates.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        Start position
    end : `~astropy.coordinates.SkyCoord`
        End position
    """

    start = ScalarSky('start')
    end = ScalarSky('end')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or {}
        self.visual = visual or {}
        self._repr_params = ('start', 'end')

    def contains(self, skycoord, wcs):
        return False

    def to_pixel(self, wcs):
        start_x, start_y = skycoord_to_pixel(self.start, wcs=wcs)
        start = PixCoord(start_x, start_y)
        end_x, end_y = skycoord_to_pixel(self.end, wcs=wcs)
        end = PixCoord(end_x, end_y)
        return LinePixelRegion(start, end)
