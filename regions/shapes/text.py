# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from .point import PointPixelRegion, PointSkyRegion
from ..core import PixCoord

__all__ = ['TextSkyRegion', 'TextPixelRegion']


class TextPixelRegion(PointPixelRegion):
    """
    A Text string in pixel coordinates

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the leftmost point of the text string
    text : str
        The text string .
    """

    def __init__(self, center, text, meta=None, visual=None):

        super(PointPixelRegion).__init__(center, meta, visual)
        self.text = text
        self._repr_params = ['Text', self.text]

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return TextSkyRegion(center, self.text)

    def as_patch(self, **kwargs):
        raise NotImplementedError


class TextSkyRegion(PointSkyRegion):
    """
    A Text string in sky coordinates

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the leftmost point of the text string
    text : str
        The text string .
    """
    def __init__(self, center, text, meta=None, visual=None):

        super(PointSkyRegion).__init__(center,  meta, visual)
        self.text = text
        self._repr_params = ['Text', self.text]

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return TextPixelRegion(center, self.text, self.meta, self.visual)
