from .core import PixelRegion, SkyRegion


class CompoundPixelRegion(PixelRegion):
    """
    Represents the logical combination of two regions in pixel coordinates.
    """

    def __init__(self, region1, operator, region2):
        self.region1 = region1
        self.region2 = region2
        self.operator = operator

    def __contains__(self, pixcoord):
        raise NotImplementedError("")

    def to_mask(self, mode='center'):
        raise NotImplementedError("")

    def to_sky(self, wcs, mode='local', tolerance=None):
        raise NotImplementedError("")

    def __repr__(self):
        return "({0} {1} {2})".format(self.region1, self.operator, self.region2)


class CompoundSkyRegion(SkyRegion):
    """
    Represents the logical combination of two regions in sky coordinates.
    """

    def __init__(self, region1, operator, region2):
        self.region1 = region1
        self.region2 = region2
        self.operator = operator

    def __contains__(self, skycoord):
        raise NotImplementedError("")

    def to_pixel(self, wcs, mode='local', tolerance=None):
        raise NotImplementedError("")

    def __repr__(self):
        return "({0} {1} {2})".format(self.region1, self.operator, self.region2)
