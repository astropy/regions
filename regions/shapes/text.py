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
    text : `str`
        The text string.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.

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

    def __init__(self, center, text, meta=None, visual=None):

        super(TextPixelRegion, self).__init__(center, meta, visual)
        self.text = text
        self._repr_params = ('text',)

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs=wcs)
        return TextSkyRegion(center, self.text)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Matplotlib Text object for this region (`matplotlib.text.Text`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        kwargs : `dict`
            All keywords that a `~matplotlib.text.Text` object accepts

        Returns
        -------
        text : `~matplotlib.text.Text`
            Matplotlib Text object.
        """
        from matplotlib.text import Text

        mpl_params = self.mpl_properties_default('text')
        mpl_params.update(kwargs)
        text = Text(self.center.x - origin[0], self.center.y - origin[1],
                    self.text, **mpl_params)

        return text


class TextSkyRegion(PointSkyRegion):
    """
    A Text string in sky coordinates

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the leftmost point of the text string
    text : `str`
        The text string.
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    def __init__(self, center, text, meta=None, visual=None):

        super(TextSkyRegion, self).__init__(center,  meta, visual)
        self.text = text
        self._repr_params = ('text',)

    def to_pixel(self, wcs):
        center_x, center_y = skycoord_to_pixel(self.center, wcs=wcs)
        center = PixCoord(center_x, center_y)
        return TextPixelRegion(center, self.text, self.meta, self.visual)
