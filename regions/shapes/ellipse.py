from astropy import units as u

from ..core import PixelRegion, SkyRegion


class EllipsePixelRegion(PixelRegion):
    """
    An ellipse in pixel coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the ellipse.
    minor : float
        The minor radius of the ellipse
    major : float
        The major radius of the ellipse
    angle : :class:`~astropy.units.Quantity`
        The rotation of the ellipse. If set to zero (the default), the major
        axis is lined up with the x axis.
    """

    def __init__(self, center, minor, major, angle=0. * u.deg, meta=None,
                 visual=None):
        # TODO: use quantity_input to check that angle is an angle
        self.center = center
        self.minor = minor
        self.major = major
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

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



class EllipseSkyRegion(SkyRegion):
    """
    An ellipse in sky coordinates.

    Parameters
    ----------
    center : :class:`~regions.core.pixcoord.PixCoord`
        The position of the center of the ellipse.
    minor : :class:`~astropy.units.Quantity`
        The minor radius of the ellipse
    major : :class:`~astropy.units.Quantity`
        The major radius of the ellipse
    angle : :class:`~astropy.units.Quantity`
        The rotation of the ellipse. If set to zero (the default), the major
        axis is lined up with the longitude axis of the celestial coordinates.
    """

    def __init__(self, center, minor, major, angle=0. * u.deg, meta=None, visual=None):
        # TODO: use quantity_input to check that height, width, and angle are angles
        self.center = center
        self.minor = minor
        self.major = major
        self.angle = angle
        self.meta = meta or {}
        self.visual = visual or {}

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
