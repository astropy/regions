from ..core import PixelRegion, SkyRegion


class PolygonPixelRegion(PixelRegion):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    vertices : :class:`~regions.core.pixcoord.PixCoord`
        The vertices of the polygon
    """

    def __init__(self, vertices):
        # TODO: test that vertices is a 1D PixCoord
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



class PolygonSkyRegion(SkyRegion):
    """
    A polygon in sky coordinates.

    Parameters
    ----------
    vertices : :class:`~regions.core.pixcoord.PixCoord`
        The vertices of the polygon
    """

    def __init__(self, vertices):
        # TODO: test that vertices is a 1D SkyCoord
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
