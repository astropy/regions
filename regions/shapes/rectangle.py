from ..core import PixelRegion, SkyRegion


class RectanglePixelRegion(PixelRegion):
    """
    An rectangle in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the rectangle.
    height : float
        The height of the rectangle
    width : float
        The width of the rectangle
    angle : :class:`~astropy.units.Quantity`
        The rotation of the rectangle. If set to zero (the default), the width
        is lined up with the x axis.
    """

    def __init__(self, vertices):
        self.vertices = vertices

    @property
    def area(self):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def __contains__(self, pixcoord):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def to_shapely(self):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def to_sky(self, wcs, mode='local', tolerance=None):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def to_mask(self, mode='center'):
        # TOOD: needs to be implemented
        raise NotImplementedError("")



class RectangleSkyRegion(SkyRegion):
    """
    An rectangle in sky coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the rectangle.
    height : :class:`~astropy.units.Quantity`
        The height radius of the rectangle
    width : :class:`~astropy.units.Quantity`
        The width radius of the rectangle
    angle : :class:`~astropy.units.Quantity`
        The rotation of the rectangle. If set to zero (the default), the width
        is lined up with the longitude axis of the celestial coordinates.
    """

    def __init__(self, vertices):
        self.vertices = vertices

    @property
    def area(self):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def __contains__(self, skycoord):
        # TOOD: needs to be implemented
        raise NotImplementedError("")

    def to_pixel(self, wcs, mode='local', tolerance=None):
        # TOOD: needs to be implemented
        raise NotImplementedError("")
