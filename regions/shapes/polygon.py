# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines polygon regions in both pixel and sky coordinates.
"""

import astropy.units as u
import numpy as np

from regions._geometry import polygonal_overlap_grid
from regions._geometry.pnpoly import points_in_polygon
from regions.core.attributes import (OneDPixCoord, OneDSkyCoord,
                                     PositiveScalar, RegionMetaDescr,
                                     RegionVisualDescr, ScalarAngle,
                                     ScalarPixCoord)
from regions.core.bounding_box import RegionBoundingBox
from regions.core.core import PixelRegion, SkyRegion
from regions.core.mask import RegionMask
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord

__all__ = ['PolygonPixelRegion', 'RegularPolygonPixelRegion',
           'PolygonSkyRegion']


class PolygonPixelRegion(PixelRegion):
    """
    A polygon in pixel coordinates.

    Parameters
    ----------
    vertices : `~regions.PixCoord`
        The vertices of the polygon.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    origin : `~regions.PixCoord`, optional
        The origin for polynomial vertices. Using this keyword allows
        ``vertices`` to be specified relative to an origin pixel
        coordinate.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, PolygonPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        x, y = [45, 45, 55, 60], [75, 70, 65, 75]
        vertices = PixCoord(x=x, y=y)
        reg = PolygonPixelRegion(vertices=vertices)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Polygon')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(30, 80)
        ax.set_ylim(50, 80)
        ax.set_aspect('equal')
    """

    _params = ('vertices',)
    _mpl_artist = 'Patch'
    vertices = OneDPixCoord('The vertices of the polygon as a |PixCoord| '
                            'array.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, vertices, meta=None, visual=None, origin=None):
        self._vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()
        if origin is None:
            origin = PixCoord(0, 0)
        self.origin = origin
        self.vertices = vertices + origin

    @property
    def area(self):
        """Return area of polygon computed by the shoelace formula."""

        # See https://stackoverflow.com/questions/24467972

        # Use offsets to improve numerical precision
        x_ = self.vertices.x - self.vertices.x.mean()
        y_ = self.vertices.y - self.vertices.y.mean()

        # Shoelace formula; for our case where the start vertex is
        # not duplicated at the end, index to avoid an array copy.
        indices = np.arange(len(x_)) - 1
        return 0.5 * abs(np.dot(x_[indices], y_) - np.dot(y_[indices], x_))

    @property
    def centroid(self):
        """Return centroid (centre of mass) of polygon."""

        # See http://paulbourke.net/geometry/polygonmesh/
        #     https://www.ma.ic.ac.uk/~rn/centroid.pdf

        # Use vertex position offsets from mean to improve numerical precision;
        # for a triangle the mean already locates the centroid.
        x0 = self.vertices.x.mean()
        y0 = self.vertices.y.mean()

        if len(self.vertices) == 3:
            return PixCoord(x0, y0)

        x_ = self.vertices.x - x0
        y_ = self.vertices.y - y0
        indices = np.arange(len(x_)) - 1

        xs = x_[indices] + x_
        ys = y_[indices] + y_
        dxy = x_[indices] * y_ - y_[indices] * x_
        scl = 1. / (6 * self.area)

        return PixCoord(np.dot(xs, dxy) * scl + x0, np.dot(ys, dxy) * scl + y0)

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, 'pixcoord')
        x = np.atleast_1d(np.asarray(pixcoord.x, dtype=float))
        y = np.atleast_1d(np.asarray(pixcoord.y, dtype=float))
        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        shape = x.shape
        mask = points_in_polygon(x.flatten(), y.flatten(), vx, vy).astype(bool)
        in_poly = mask.reshape(shape)
        if self.meta.get('include', True):
            return in_poly
        else:
            return np.logical_not(in_poly)

    def to_sky(self, wcs):
        vertices_sky = wcs.pixel_to_world(self.vertices.x, self.vertices.y)
        return PolygonSkyRegion(vertices=vertices_sky, meta=self.meta.copy(),
                                visual=self.visual.copy())

    @property
    def bounding_box(self):
        xmin = self.vertices.x.min()
        xmax = self.vertices.x.max()
        ymin = self.vertices.y.min()
        ymax = self.vertices.y.max()
        return RegionBoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):
        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        use_exact = 0 if mode == 'subpixels' else 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at origin
        xmin = float(bbox.ixmin) - 0.5
        xmax = float(bbox.ixmax) - 0.5
        ymin = float(bbox.iymin) - 0.5
        ymax = float(bbox.iymax) - 0.5

        vx = np.asarray(self.vertices.x, dtype=float)
        vy = np.asarray(self.vertices.y, dtype=float)

        fraction = polygonal_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                          vx, vy, use_exact, subpixels)

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Polygon`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Polygon`. These keywords will override
            any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Polygon`
            A matplotlib polygon patch.
        """
        from matplotlib.patches import Polygon

        xy = np.vstack([self.vertices.x - origin[0],
                        self.vertices.y - origin[1]]).transpose()

        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Polygon(xy=xy, **mpl_kwargs)

    def _update_from_mpl_selector(self, verts, *args, **kwargs):
        """Set position and orientation from selector properties."""
        # Polygon selector calls ``callback(self.verts)``.

        self.vertices = PixCoord(*np.array(verts).T)

        if getattr(self, '_mpl_selector_callback', None) is not None:
            self._mpl_selector_callback(self)

    def as_mpl_selector(self, ax, active=True, sync=True, callback=None, **kwargs):
        """
        A matplotlib editable widget for this region
        (`matplotlib.widgets.PolygonSelector`).

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
            Additional keyword arguments that are passed to
            `matplotlib.widgets.PolygonSelector`.

        Returns
        -------
        selector : `matplotlib.widgets.PolygonSelector`
            The matplotlib selector.

        Notes
        -----
        Once a selector has been created, you will need to keep a
        reference to it until you no longer need it. In addition,
        you can enable/disable the selector at any point by calling
        ``selector.set_active(True)`` or ``selector.set_active(False)``.
        """
        from matplotlib.widgets import PolygonSelector
        import matplotlib._version
        _mpl_version = getattr(matplotlib._version, 'version', None)
        if _mpl_version is None:
            _mpl_version = matplotlib._version.get_versions()['version']

        if hasattr(self, '_mpl_selector'):
            raise Exception('Cannot attach more than one selector to a region.')

        if not hasattr(PolygonSelector, '_scale_polygon'):
            raise NotImplementedError('Rescalable ``PolygonSelector`` widgets are not '
                                      f'yet supported with matplotlib {_mpl_version}.')

        if sync:
            sync_callback = self._update_from_mpl_selector
        else:
            def sync_callback(*args, **kwargs):
                pass

        self._mpl_selector = PolygonSelector(
            ax, sync_callback, draw_bounding_box=True,
            props={'color': self.visual.get('color', 'black'),
                   'linewidth': self.visual.get('linewidth', 1),
                   'linestyle': self.visual.get('linestyle', 'solid')})

        self._mpl_selector.verts = list(zip(self.vertices.x, self.vertices.y))
        self._mpl_selector.set_active(active)
        self._mpl_selector_callback = callback

        if sync and self._mpl_selector_callback is not None:
            self._mpl_selector_callback(self)

        return self._mpl_selector

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
        region : `PolygonPixelRegion`
            The rotated region (which is an independent copy).
        """
        vertices = self.vertices.rotate(center, angle)
        return self.copy(vertices=vertices)


class RegularPolygonPixelRegion(PolygonPixelRegion):
    """
    A regular polygon in pixel coordinates.

    .. note::

        This class will be serialized as a generic polygon
        region, thus when read back in it will produce a
        `~regions.PolygonPixelRegion` object instead of a
        `~regions.RegularPolygonPixelRegion` object.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the polygon.
    nvertices : `~regions.PixCoord`
        The number of polygon vertices (or sides).
    radius : float
        The distance from the center to any vertex. This is also known
        as the circumradius.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the polygon, measured anti-clockwise. If
        set to zero (the default), the polygon will point "up" following
        the `matplotlib.patches.RegularPolygon` convention.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Attributes
    ----------
    side_length : float
        The side length.
    inradius : float
        The radius of the largest circle contained entirely within the
        polygon. This value is identical to the length of a line segment
        from the polygon center to the midpoint of one of its sides
        (known as as the apothem).
    perimeter : float
        The polygon perimeter.
    interior_angle : float
        The polygon interior angle, which is the angle at each vertex on
        the inside of the polygon.
    exterior_angle : float
        The polygon exterior angle, which is an angle at each vertex on
        the outside of the polygon.

    Examples
    --------
    .. plot::
        :include-source:

        import astropy.units as u
        import matplotlib.pyplot as plt
        from regions import PixCoord, RegularPolygonPixelRegion

        fig, ax = plt.subplots(1, 1)
        center = PixCoord(x=50, y=50)
        reg1 = RegularPolygonPixelRegion(center, 6, 15)
        reg1.plot(edgecolor='red', lw=2)

        center = PixCoord(x=25, y=25)
        reg2 = RegularPolygonPixelRegion(center, 3, 15)
        reg2.plot(edgecolor='green', lw=2)

        center = PixCoord(x=25, y=75)
        reg3 = RegularPolygonPixelRegion(center, 3, 15, angle=25*u.deg)
        reg3.plot(edgecolor='orange', lw=2)

        center = PixCoord(x=75, y=75)
        reg4 = RegularPolygonPixelRegion(center, 8, 15)
        reg4.plot(edgecolor='blue', lw=2)

        center = PixCoord(x=75, y=25)
        reg5 = RegularPolygonPixelRegion(center, 5, 15)
        reg5.plot(edgecolor='magenta', lw=2)

        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_aspect('equal')
    """

    _params = ('center', 'nvertices', 'radius', 'angle')
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    nvertices = PositiveScalar('The number of polygon vertices.')
    radius = PositiveScalar('The distance from the center to any vertex in '
                            'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, nvertices, radius, angle=0. * u.deg,
                 meta=None, visual=None):

        if nvertices < 3:
            raise ValueError('nvertices must be >= 3')
        self.center = center
        self.nvertices = nvertices
        self.radius = radius
        self.angle = angle

        super().__init__(self._calc_vertices(), meta=meta, visual=visual)

        self.side_length = 2. * self.radius * np.sin(np.pi / self.nvertices)
        self.inradius = self.radius * np.cos(np.pi / self.nvertices)
        self.perimeter = self.side_length * self.nvertices
        self.interior_angle = ((self.nvertices - 2) / self.nvertices
                               * 180) << u.deg
        self.exterior_angle = 360. / self.nvertices << u.deg

    def _calc_vertices(self):
        # uses the matplotlib convention that the polygon always points "up"
        theta = ((2. * np.pi / self.nvertices * np.arange(self.nvertices)
                 + (np.pi / 2)) << u.radian) + self.angle
        xvert = self.radius * np.cos(theta)
        yvert = self.radius * np.sin(theta)
        return self.center + PixCoord(xvert, yvert)

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
        region : `PolygonPixelRegion`
            The rotated region (which is an independent copy).
        """
        center = self.center.rotate(center, angle)
        angle = self.angle + angle
        return self.copy(center=center, angle=angle)

    def to_polygon(self):
        """
        Return a `PolygonPixelRegion` of this region.
        """
        return PolygonPixelRegion(self.vertices, origin=PixCoord(0, 0),
                                  meta=self.meta.copy(),
                                  visual=self.visual.copy())


class PolygonSkyRegion(SkyRegion):
    """
    A polygon defined using vertices in sky coordinates.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('vertices',)
    vertices = OneDSkyCoord('The vertices of the polygon as a |SkyCoord| '
                            'array.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, vertices, meta=None, visual=None):
        self.vertices = vertices
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        x, y = wcs.world_to_pixel(self.vertices)
        vertices_pix = PixCoord(x, y)
        return PolygonPixelRegion(vertices_pix, meta=self.meta.copy(),
                                  visual=self.visual.copy())
