# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines rectangular regions in both pixel and sky coordinates.
"""

from astropy.coordinates import Angle
import astropy.units as u
from astropy.wcs.utils import pixel_to_skycoord
import numpy as np

from ..core.attributes import (ScalarPix, ScalarLength, QuantityLength,
                               ScalarSky)
from ..core.bounding_box import BoundingBox
from ..core.core import PixelRegion, SkyRegion
from ..core.mask import RegionMask
from ..core.metadata import RegionMeta, RegionVisual
from ..core.pixcoord import PixCoord
from .._geometry import rectangular_overlap_grid
from .._utils.wcs_helpers import pixel_scale_angle_at_skycoord

from .polygon import PolygonPixelRegion

__all__ = ['RectanglePixelRegion', 'RectangleSkyRegion']


class RectanglePixelRegion(PixelRegion):
    """
    A rectangle in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the rectangle.
    width : float
        The width of the rectangle (before rotation) in pixels.
    height : float
        The height of the rectangle (before rotation) in pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangle, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the x
        axis.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from astropy.coordinates import Angle
        from regions import PixCoord, RectanglePixelRegion
        import matplotlib.pyplot as plt

        x, y = 15, 10
        width, height = 8, 5
        angle = Angle(30, 'deg')

        fig, ax = plt.subplots(1, 1)

        center = PixCoord(x=x, y=y)
        reg = RectanglePixelRegion(center=center, width=width,
                                   height=height, angle=angle)
        patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
        ax.add_patch(patch)

        plt.xlim(0, 30)
        plt.ylim(0, 20)
        ax.set_aspect('equal')
    """

    _params = ('center', 'width', 'height', 'angle')
    center = ScalarPix('center')
    width = ScalarLength('width')
    height = ScalarLength('height')
    angle = QuantityLength('angle')
    mpl_artist = 'Patch'

    def __init__(self, center, width, height, angle=0 * u.deg, meta=None,
                 visual=None):
        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
        return self.width * self.height

    def contains(self, pixcoord):
        cos_angle = np.cos(self.angle)
        sin_angle = np.sin(self.angle)
        dx = pixcoord.x - self.center.x
        dy = pixcoord.y - self.center.y
        dx_rot = cos_angle * dx + sin_angle * dy
        dy_rot = sin_angle * dx - cos_angle * dy
        in_rect = ((np.abs(dx_rot) < self.width * 0.5)
                   & (np.abs(dy_rot) < self.height * 0.5))
        if self.meta.get('include', True):
            return in_rect
        else:
            return np.logical_not(in_rect)

    def to_sky(self, wcs):
        center = pixel_to_skycoord(self.center.x, self.center.y, wcs)
        _, pixscale, north_angle = pixel_scale_angle_at_skycoord(center, wcs)
        width = Angle(self.width * u.pix * pixscale, 'arcsec')
        height = Angle(self.height * u.pix * pixscale, 'arcsec')
        angle = self.angle - (north_angle - 90 * u.deg)
        return RectangleSkyRegion(center, width, height, angle=angle,
                                  meta=self.meta, visual=self.visual)

    @property
    def bounding_box(self):
        w2 = self.width / 2.
        h2 = self.height / 2.
        cos_angle = np.cos(self.angle)  # self.angle is a Quantity
        sin_angle = np.sin(self.angle)  # self.angle is a Quantity
        dx1 = abs(w2 * cos_angle - h2 * sin_angle)
        dy1 = abs(w2 * sin_angle + h2 * cos_angle)
        dx2 = abs(w2 * cos_angle + h2 * sin_angle)
        dy2 = abs(w2 * sin_angle - h2 * cos_angle)
        dx = max(dx1, dx2)
        dy = max(dy1, dy2)

        xmin = self.center.x - dx
        xmax = self.center.x + dx
        ymin = self.center.y - dy
        ymax = self.center.y + dy

        return BoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):
        # NOTE: assumes this class represents a single rectangle

        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at
        # origin
        xmin = float(bbox.ixmin) - 0.5 - self.center.x
        xmax = float(bbox.ixmax) - 0.5 - self.center.x
        ymin = float(bbox.iymin) - 0.5 - self.center.y
        ymax = float(bbox.iymax) - 0.5 - self.center.y

        if mode == 'subpixels':
            use_exact = 0
        else:
            use_exact = 1

        fraction = rectangular_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                            self.width, self.height,
                                            self.angle.to(u.rad).value,
                                            use_exact, subpixels)

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Rectangle`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Rectangle`. These keywords will
            override any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Rectangle`
            A matplotlib rectangle patch.
        """
        from matplotlib.patches import Rectangle

        xy = self._lower_left_xy()
        xy = xy[0] - origin[0], xy[1] - origin[1]
        width = self.width
        height = self.height
        # matplotlib expects rotation in degrees (anti-clockwise)
        angle = self.angle.to('deg').value

        mpl_kwargs = self.visual.define_mpl_kwargs(self.mpl_artist)
        mpl_kwargs.update(kwargs)

        return Rectangle(xy=xy, width=width, height=height,
                         angle=angle, **mpl_kwargs)

    def _update_from_mpl_selector(self, *args, **kwargs):
        xmin, xmax, ymin, ymax = self._mpl_selector.extents
        self.center = PixCoord(x=0.5 * (xmin + xmax),
                               y=0.5 * (ymin + ymax))
        self.width = (xmax - xmin)
        self.height = (ymax - ymin)
        self.angle = 0. * u.deg
        if self._mpl_selector_callback is not None:
            self._mpl_selector_callback(self)

    def as_mpl_selector(self, ax, active=True, sync=True, callback=None,
                        **kwargs):
        """
        A matplotlib editable widget for this region
        (`matplotlib.widgets.RectangleSelector`).

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            The matplotlib axes to add the selector to.
        active : bool, optional
            Whether the selector should be active by default.
        sync : bool, optional
            If `True` (the default), the region will be kept in
            sync with the selector. Otherwise, the selector will be
            initialized with the values from the region but the two will
            then be disconnected.
        callback : callable, optional
            If specified, this function will be called every time the
            region is updated. This only has an effect if ``sync`` is
            `True`. If a callback is set, it is called for the first
            time once the selector has been created.
        **kwargs : dict
            Additional keyword arguments are passed to
            `matplotlib.widgets.RectangleSelector`.

        Returns
        -------
        selector : `matplotlib.widgets.RectangleSelector`
            The matplotlib selector.

        Notes
        -----
        Once a selector has been created, you will need to keep a
        reference to it until you no longer need it. In addition,
        you can enable/disable the selector at any point by calling
        ``selector.set_active(True)`` or ``selector.set_active(False)``.
        """
        from matplotlib.widgets import RectangleSelector

        if hasattr(self, '_mpl_selector'):
            raise Exception('Cannot attach more than one selector to a '
                            'region.')

        if self.angle.value != 0:
            raise NotImplementedError('Cannot create matplotlib selector for '
                                      'rotated rectangle.')

        if sync:
            sync_callback = self._update_from_mpl_selector
        else:
            def sync_callback(*args, **kwargs):
                pass

        rectprops = kwargs.pop('rectprops', {'edgecolor': self.visual.get('color', 'black'),
                                             'facecolor': 'none',
                                             'linewidth': self.visual.get('linewidth', 1),
                                             'linestyle': self.visual.get('linestyle', 'solid')})

        self._mpl_selector = RectangleSelector(ax, sync_callback, interactive=True,
                                               rectprops=rectprops, **kwargs)

        self._mpl_selector.extents = (self.center.x - self.width / 2,
                                      self.center.x + self.width / 2,
                                      self.center.y - self.height / 2,
                                      self.center.y + self.height / 2)
        self._mpl_selector.set_active(active)
        self._mpl_selector_callback = callback

        if sync and self._mpl_selector_callback is not None:
            self._mpl_selector_callback(self)

        return self._mpl_selector

    @property
    def corners(self):
        """
        Return the x, y coordinate pairs that define the corners.
        """
        corners = [(-self.width / 2, -self.height / 2),
                   (self.width / 2, -self.height / 2),
                   (self.width / 2, self.height / 2),
                   (-self.width / 2, self.height / 2),
                   ]
        rotmat = [[np.cos(self.angle), np.sin(self.angle)],
                  [-np.sin(self.angle), np.cos(self.angle)]]

        return np.dot(corners, rotmat) + np.array([self.center.x,
                                                   self.center.y])

    def to_polygon(self):
        """
        Return a 4-sided polygon equivalent to this rectangle.
        """
        x, y = self.corners.T
        vertices = PixCoord(x=x, y=y)
        return PolygonPixelRegion(vertices=vertices, meta=self.meta,
                                  visual=self.visual)

    def _lower_left_xy(self):
        """
        Compute lower left ``xy`` pixel position.

        This is used for the conversion to matplotlib in ``as_artist``.

        Taken from
        http://photutils.readthedocs.io/en/latest/_modules/photutils/ape
        rture/rectangle.html#RectangularAperture.plot.
        """
        hw = self.width / 2.
        hh = self.height / 2.
        sint = np.sin(self.angle)
        cost = np.cos(self.angle)
        dx = (hh * sint) - (hw * cost)
        dy = -(hh * cost) - (hw * sint)
        x = self.center.x + dx
        y = self.center.y + dy
        return x, y

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `RectanglePixelRegion`
            The rotated region (which is an independent copy).
        """
        center = self.center.rotate(center, angle)
        angle = self.angle + angle
        return self.copy(center=center, angle=angle)


class RectangleSkyRegion(SkyRegion):
    """
    A rectangle in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangle.
    width : `~astropy.units.Quantity`
        The width of the rectangle (before rotation) as an angle.
    height : `~astropy.units.Quantity`
        The height of the rectangle (before rotation) as an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangle, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the
        longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta`, optional
        A dictionary that stores the meta attributes of this region.
    visual : `~regions.RegionVisual`, optional
        A dictionary that stores the visual meta attributes of this
        region.
    """

    _params = ('center', 'width', 'height', 'angle')
    center = ScalarSky('center')
    width = QuantityLength('width')
    height = QuantityLength('height')
    angle = QuantityLength('angle')

    def __init__(self, center, width, height, angle=0 * u.deg, meta=None,
                 visual=None):
        self.center = center
        self.width = width
        self.height = height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        center, pixscale, north_angle = pixel_scale_angle_at_skycoord(
            self.center, wcs)
        width = (self.width / pixscale).to(u.pix).value
        height = (self.height / pixscale).to(u.pix).value
        angle = self.angle + (north_angle - 90 * u.deg)
        return RectanglePixelRegion(center, width, height, angle=angle,
                                    meta=self.meta, visual=self.visual)
