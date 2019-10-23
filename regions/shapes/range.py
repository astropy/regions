# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.coordinates import Angle, SkyCoord
from regions.core.attributes import QuantityLength
from ..core import SkyRegion, PixelRegion, PixCoord

__all__ = ["RangePixelRegion", "RangeSphericalRegion"]

inf = Angle(float("inf"), "deg")


# TODO: needed?
class RangePixelRegion(PixelRegion):
    """
    A range in pixel coordinates.

    See `RangeSphericalRegion`.
    """


class RangeSphericalRegion(SkyRegion):
    """
    A range in longitude and latitude.

    Use case: http://www.ivoa.net/documents/SIA/20151223/REC-SIA-2.0-20151223.html#toc12

    Parameters
    ----------
    lon_min, lon_max : `~astropy.coordinates.Angle`
        Range in longitude
    lat_min, lat_max : `~astropy.coordinates.Angle`
        Range in latitude
    frame : `~astropy.coordinates.BaseCoordinateFrame`
        Sky frame
    meta : `~regions.RegionMeta` object, optional
        A dictionary which stores the meta attributes of this region.
    visual : `~regions.RegionVisual` object, optional
        A dictionary which stores the visual meta attributes of this region.
    """
    # TODO: should we use these properties, or SkyCoord for min / max,
    # or a center-based representation?
    _params = ('lon_min', 'lon_max', 'lat_min', 'lat_max')
    lon_min = QuantityLength('lon_min')
    lon_max = QuantityLength('lon_max')
    lat_min = QuantityLength('lat_min')
    lat_max = QuantityLength('lat_max')

    def __init__(self, lon_min=-inf, lon_max=inf, lat_min=-inf, lat_max=inf, frame="icrs", meta=None, visual=None):
        self.lon_min = lon_min
        self.lon_max = lon_max
        self.lat_min = lat_min
        self.lat_max = lat_max
        # TODO: Is there a helper function to make a frame without a
        self.frame = SkyCoord(0, 0, unit="deg", frame=frame).frame
        self.meta = meta or {}
        self.visual = visual or {}

    def to_pixel(self, wcs):
        sky_min = SkyCoord(self.lon_min, self.lat_min, frame=self.frame)
        sky_max = SkyCoord(self.lon_max, self.lat_max, frame=self.frame)

        pix_min = PixCoord.from_sky(sky_min, self.wcs)
        pix_max = PixCoord.from_sky(sky_max, self.wcs)
        return RangePixelRegion(pix_min.x, pix_max.x, pix_min.y, pix_max.y, self.meta, self.visual)

    def contains(self, skycoord, wcs=None):
        skycoord = skycoord.transform_to(self.frame)
        return (
            (self.lon_min <= skycoord.data.lon) &
            (self.lon_max >= skycoord.data.lon) &
            (self.lat_min <= skycoord.data.lat) &
            (self.lat_max >= skycoord.data.lat)
        )
