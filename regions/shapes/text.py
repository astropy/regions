# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines text regions in both pixel and sky coordinates.
"""

from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

from .point import PointPixelRegion, PointSkyRegion
from ..core import PixCoord

__all__ = ['TextSkyRegion', 'TextPixelRegion']


class TextPixelRegion(PointPixelRegion):
    """
    A text string in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the leftmost point of the text string.
    text : str
        The text string.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, TextPixelRegion, RegionVisual
        import matplotlib.pyplot as plt

        x, y = 15, 10
        textangle = '30'

        fig, ax = plt.subplots(1, 1)

        center = PixCoord(x=x, y=y)
        reg = TextPixelRegion(center=center, text="Hello World!",
                              visual=RegionVisual(textangle=textangle))
        reg.plot(ax=ax)

        plt.xlim(10, 30)
        plt.ylim(2.5, 20)
        ax.set_aspect('equal')
    """

    _params = ('center', 'text')
    mpl_artist = 'Text'

    def __init__(self, center, text, meta=None, visual=None):
        super().__init__(center, meta, visual)
        self.text = text

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return TextSkyRegion(center, self.text, meta=self.meta.copy(),
                             visual=self.visual.copy())

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib Text object for this region
        (`matplotlib.text.Text`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by `~matplotlib.text.Text`.
            These keywords will override any visual meta attributes of
            this region.

        Returns
        -------
        artist : `~matplotlib.text.Text`
            A matplotlib Text object.
        """
        from matplotlib.text import Text

        mpl_kwargs = self.visual.define_mpl_kwargs(self.mpl_artist)
        mpl_kwargs.update(kwargs)

        return Text(self.center.x - origin[0], self.center.y - origin[1],
                    self.text, **mpl_kwargs)


class TextSkyRegion(PointSkyRegion):
    """
    A text string in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the leftmost point of the text string.
    text : str
        The text string.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('center', 'text')

    def __init__(self, center, text, meta=None, visual=None):
        super().__init__(center, meta, visual)
        self.text = text

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return TextPixelRegion(center, self.text, meta=self.meta.copy(),
                               visual=self.visual.copy())
