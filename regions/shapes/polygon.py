# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines polygon regions in both pixel and sky coordinates.
"""
import operator

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.stats import circmean

from regions._geometry import polygonal_overlap_grid
from regions._geometry.pnpoly import points_in_polygon
from regions._utils.spherical_helpers import (
    cross_product_skycoord2skycoord, cross_product_sum_skycoord2skycoord,
    discretize_all_edge_boundaries, get_edge_raw_lonlat_bounds_circ_edges)
from regions.core.attributes import (OneDPixCoord, OneDSkyCoord,
                                     PositiveScalar, RegionMetaDescr,
                                     RegionVisualDescr, ScalarAngle,
                                     ScalarPixCoord)
from regions.core.bounding_box import RegionBoundingBox
from regions.core.compound import CompoundSphericalSkyRegion
from regions.core.core import PixelRegion, SkyRegion, SphericalSkyRegion
from regions.core.mask import RegionMask
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord

__all__ = ['PolygonPixelRegion', 'RegularPolygonPixelRegion',
           'PolygonSkyRegion', 'PolygonSphericalSkyRegion']


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
        # See https://stackoverflow.com/questions/24467972

        # Use offsets to improve numerical precision
        x_ = self.vertices.x - self.vertices.x.mean()
        y_ = self.vertices.y - self.vertices.y.mean()

        # Shoelace formula, for our case where the start vertex
        # isn't duplicated at the end, written to avoid an array copy
        area_main = np.dot(x_[:-1], y_[1:]) - np.dot(y_[:-1], x_[1:])
        area_last = x_[-1] * y_[0] - y_[-1] * x_[0]
        return 0.5 * np.abs(area_main + area_last)

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

    def to_spherical_sky(self, wcs=None, include_boundary_distortions=False,
                         discretize_kwargs=None):
        if discretize_kwargs is None:
            discretize_kwargs = {}

        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires cylindrical to spherical projection (using WCS) and discretization
            # Will require implementing discretization in pixel space
            # to get correct handling of distortions.
            raise NotImplementedError

            # ### Potential solution:
            # # Leverage polygon class to_spherical_sky() functionality without
            # # distortions, as the distortions were already computed in creating
            # # that polygon approximation
            # return self.to_pixel(wcs).discretize_boundary(**discretize_kwargs).to_spherical_sky(
            #     wcs=wcs, include_boundary_distortions=False
            # )

        # TODO: ensure vertices are in CW order,
        # as implicitly a planar -> spherical polygon will not be the
        # "large" complement polygon on the sphere?

        return self.to_sky(wcs).to_spherical_sky()

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

        # Find position of pixel edges and recenter so that circle is at
        # origin
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

    def to_spherical_sky(self, wcs=None, include_boundary_distortions=False,
                         discretize_kwargs=None):
        if discretize_kwargs is None:
            discretize_kwargs = {}

        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires cylindrical to spherical projection (using WCS) and discretization
            # Will require implementing discretization in pixel space
            # to get correct handling of distortions.
            raise NotImplementedError

            # ### Potential solution:
            # # Leverage polygon class to_spherical_sky() functionality without
            # # distortions, as the distortions were already computed in creating
            # # that polygon approximation
            # return self.to_pixel(wcs).discretize_boundary(**discretize_kwargs).to_spherical_sky(
            #     wcs=wcs, include_boundary_distortions=False
            # )

        # TODO: ensure vertices are in CW order,
        # as implicitly a planar -> spherical polygon will not be the
        # "large" complement polygon on the sphere?

        return PolygonSphericalSkyRegion(
            self.vertices,
            self.meta.copy(),
            self.visual.copy()
        )


class PolygonSphericalSkyRegion(SphericalSkyRegion):
    """
    A spherical polygon defined using vertices in sky coordinates.

    Internal "contains" logic could be changed to
    leverage spherical-geometry package functionality.

    This region is very much akin to the `~regions.PolygonPixelRegion`
    class (and borrows internal attributes following the same structure).

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the polygon. Assumed to be in CW order.
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

    @property
    def _edge_circs(self):
        """
        Get list of the great circles defining the polygon boundaries.
        """
        from .circle import CircleSphericalSkyRegion
        gcs = []
        for i in range(len(self.vertices)):
            # verts are in CW order: Cross product to get bounding great circle centers
            # Compute GCs and stack into a compound set:
            c_gc = cross_product_skycoord2skycoord(self.vertices[i - 1], self.vertices[i])
            gcs.append(CircleSphericalSkyRegion(c_gc, 90 * u.deg))
        return gcs

    @property
    def _compound_region(self):
        # Need N great circles to define boundaries for an N-sided polygon.
        # verts are in CW order: Cross product to get bounding great circle centers
        # Compute GCs and stack into a compound set:
        compreg = None
        gcs = self._edge_circs
        for gc in gcs:
            if compreg is None:
                compreg = gc
            else:
                compreg = CompoundSphericalSkyRegion(
                    compreg, gc, operator.and_, self.meta, self.visual,
                )

        return compreg

    @property
    def centroid(self):
        """
        Region centroid.

        Defined as the point equidistant from all vertices.

        However, if this is not contained within the polygon, instead
        use the average of the vertices' positions.
        """
        # Calculate from cross products of vertices with cartesian representation
        # verts are in CW order:
        # Minimum distance
        centroid_mindist = self.centroid_mindist

        if not self.contains(centroid_mindist):
            return self.centroid_avg

        return centroid_mindist

    @property
    def centroid_mindist(self):
        """
        Region centroid, defined as the point equidistant from all
        vertices (minimum distance).
        """
        # Calculate from cross products of vertices with cartesian representation
        # verts are in CW order:
        # Minimum distance
        return cross_product_sum_skycoord2skycoord(self.vertices)

    @property
    def centroid_avg(self):
        """
        Region centroid, taking the average of the vertices'
        coordinates.
        """
        verts_sph = self.vertices.represent_as('spherical')

        lons = verts_sph.lon
        lats = verts_sph.lat

        lon = circmean(lons)
        lat = np.average(lats)
        return SkyCoord(lon, lat, frame=self.frame)

    @property
    def bounding_circle(self):
        from .circle import CircleSphericalSkyRegion
        cent = self.centroid
        seps = cent.separation(self.vertices)
        return CircleSphericalSkyRegion(center=cent, radius=np.max(seps))

    @property
    def bounding_lonlat(self):
        lons_arr, lats_arr = get_edge_raw_lonlat_bounds_circ_edges(
            self.vertices, self.centroid, self._edge_circs
        )

        lons_arr, lats_arr = self._validate_lonlat_bounds(lons_arr, lats_arr)

        return lons_arr, lats_arr

    def contains(self, coord):
        return self._compound_region.contains(coord)

    def transform_to(self, frame, merge_attributes=True):
        frame = self._validate_frame(frame)

        # Only center transforms, radii preserved
        verts_transf = self.vertices.transform_to(frame, merge_attributes=merge_attributes)

        return PolygonSphericalSkyRegion(
            verts_transf,
            self.meta.copy(),
            self.visual.copy()
        )

    def discretize_boundary(self, n_points=10):
        bound_verts = discretize_all_edge_boundaries(
            self.vertices, self._edge_circs, self.centroid, n_points
        )
        return PolygonSphericalSkyRegion(bound_verts)

    def to_sky(
            self,
            wcs=None,
            include_boundary_distortions=False,
            discretize_kwargs=None
    ):

        if discretize_kwargs is None:
            discretize_kwargs = {}

        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires spherical to cylindrical projection (from WCS) and discretization
            # Use to_pixel(), then apply "small angle approx" to get planar sky.
            return self.to_pixel(
                include_boundary_distortions=include_boundary_distortions,
                wcs=wcs,
                discretize_kwargs=discretize_kwargs,
            ).to_sky(wcs)

        return PolygonSkyRegion(
            self.vertices,
            meta=self.meta,
            visual=self.visual
        )

    def to_pixel(
            self,
            wcs=None,
            include_boundary_distortions=False,
            discretize_kwargs=None,
    ):

        if discretize_kwargs is None:
            discretize_kwargs = {}
        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires spherical to cylindrical projection (from WCS) and discretization
            verts = wcs.world_to_pixel(
                self.discretize_boundary(**discretize_kwargs).vertices
            )

            return PolygonPixelRegion(
                PixCoord(*verts), meta=self.meta.copy(),
                visual=self.visual.copy()
            )

        return self.to_sky().to_pixel(wcs)
