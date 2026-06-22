# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines annulus regions in both pixel and sky coordinates.
"""

import abc
import math
import operator

import astropy.units as u
import numpy as np

from regions._utils.wcs_helpers import (pixel_shape_to_sky_svd,
                                        pixel_to_sky_mean_scale,
                                        pixel_to_sky_svd_scales,
                                        sky_shape_to_pixel_svd,
                                        sky_to_pixel_mean_scale,
                                        sky_to_pixel_svd_scales)
from regions.core.attributes import (PositiveScalar, PositiveScalarAngle,
                                     RegionMetaDescr, RegionVisualDescr,
                                     ScalarAngle, ScalarPixCoord,
                                     ScalarSkyCoord)
from regions.core.compound import (CompoundPixelRegion,
                                   CompoundSphericalSkyRegion)
from regions.core.core import PixelRegion, SkyRegion, SphericalSkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord
from regions.shapes.circle import CirclePixelRegion, CircleSphericalSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion

__all__ = ['AnnulusPixelRegion', 'AnnulusSphericalSkyRegion',
           'AsymmetricAnnulusPixelRegion',
           'AsymmetricAnnulusSkyRegion',
           'CircleAnnulusPixelRegion', 'CircleAnnulusSkyRegion',
           'CircleAnnulusSphericalSkyRegion',
           'EllipseAnnulusPixelRegion', 'EllipseAnnulusSkyRegion',
           'RectangleAnnulusPixelRegion', 'RectangleAnnulusSkyRegion']


class AnnulusPixelRegion(PixelRegion, abc.ABC):
    """
    A base class for annulus pixel regions.
    """

    @property
    def _compound_region(self):
        return CompoundPixelRegion(self._inner_region, self._outer_region,
                                   operator.xor, self.meta, self.visual)

    @property
    def area(self):
        return self._outer_region.area - self._inner_region.area

    @property
    def bounding_box(self):
        return self._outer_region.bounding_box

    def contains(self, pixcoord):
        return self._compound_region.contains(pixcoord)

    def covers(self, pixcoord):
        outer_cov = self._outer_region.covers(pixcoord)
        inner_con = self._inner_region.contains(pixcoord)
        return np.logical_and(outer_cov, np.logical_not(inner_con))

    def as_artist(self, origin=(0, 0), **kwargs):
        return self._compound_region.as_artist(origin, **kwargs)

    def to_mask(self, mode='center', subpixels=5):
        return self._compound_region.to_mask(mode, subpixels)

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
        region : `~regions.AnnulusPixelRegion`
            The rotated region (which is an independent copy).
        """
        changes = {}
        changes['center'] = self.center.rotate(center, angle)
        if hasattr(self, 'angle'):
            changes['angle'] = self.angle + angle
        return self.copy(**changes)


class AnnulusSphericalSkyRegion(SphericalSkyRegion, abc.ABC):
    """
    A base class for spherical sky annulus regions.
    """

    @property
    def _compound_region(self):
        return CompoundSphericalSkyRegion(
            self._inner_region, self._outer_region,
            operator.xor, self.meta, self.visual)

    @property
    def frame(self):
        return self._inner_region.frame

    @property
    def bounding_circle(self):
        return self._outer_region.bounding_circle

    @property
    def bounding_lonlat(self):
        # Bounding lonlat of inner, outer regions
        # Check of inner region goes over the pole so
        # a more limited lat range is needed
        lons_arr, lats_arr = self._outer_region.bounding_lonlat

        # Check if shape covers either pole & modify lats arr accordingly,
        # accounting for annular geometry:
        lons_arr, lats_arr = self._validate_lonlat_bounds(
            lons_arr, lats_arr, inner_region=self._inner_region,
        )

        return lons_arr, lats_arr

    def contains(self, coord):
        return self._compound_region.contains(coord)

    def covers(self, coord):
        return self._compound_region.covers(coord)


class CircleAnnulusPixelRegion(AnnulusPixelRegion):
    """
    A circular annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    inner_radius : float
        The inner radius of the annulus in pixels.
    outer_radius : float
        The outer radius of the annulus in pixels.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from regions import CircleAnnulusPixelRegion, PixCoord

        fig, ax = plt.subplots()

        reg = CircleAnnulusPixelRegion(PixCoord(x=6, y=6),
                                       inner_radius=5.5,
                                       outer_radius=8.0)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Circle Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    _component_class = CirclePixelRegion
    _params = ('center', 'inner_radius', 'outer_radius')
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_radius = PositiveScalar('The inner radius in pixels as a float.')
    outer_radius = PositiveScalar('The outer radius in pixels as a float.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_radius, outer_radius, meta=None,
                 visual=None):
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_radius >= outer_radius:
            raise ValueError('outer_radius must be greater than inner_radius')

    @property
    def _inner_region(self):
        return self._component_class(self.center, self.inner_radius,
                                     self.meta, self.visual)

    @property
    def _outer_region(self):
        return self._component_class(self.center, self.outer_radius,
                                     self.meta, self.visual)

    def to_sky(self, wcs, *, as_ellipse=False):
        """
        Return a sky region from this pixel region.

        Parameters
        ----------
        wcs : WCS object
            A world coordinate system (WCS) transformation that
            supports the `astropy shared interface for WCS
            <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_
            (e.g., `astropy.wcs.WCS`).

        as_ellipse : bool, optional
            If `True`, return an `~regions.EllipseAnnulusSkyRegion`
            instead of a `~regions.CircleAnnulusSkyRegion`. An ellipse
            annulus is generally a better approximation when the WCS has
            distortions or different pixel scales along different axes.
            Default is `False`.

        Returns
        -------
        region : `~regions.CircleAnnulusSkyRegion` or `~regions.EllipseAnnulusSkyRegion`
            The sky region. An ellipse annulus is returned if
            ``as_ellipse`` is `True`.
        """
        if as_ellipse:
            center, scale_major, scale_minor, angle = pixel_to_sky_svd_scales(
                (self.center.x, self.center.y), wcs)
            inner_width = 2 * self.inner_radius * scale_major * u.arcsec
            outer_width = 2 * self.outer_radius * scale_major * u.arcsec
            inner_height = 2 * self.inner_radius * scale_minor * u.arcsec
            outer_height = 2 * self.outer_radius * scale_minor * u.arcsec
            # The helper returns a position angle (PA) from North;
            # regions measures the angle from the RA axis (90 deg
            # offset).
            angle = (angle + 90 * u.deg).wrap_at(360 * u.deg)
            return EllipseAnnulusSkyRegion(
                center, inner_width, outer_width,
                inner_height, outer_height, angle=angle,
                meta=self.meta.copy(), visual=self.visual.copy())

        center, mean_scale = pixel_to_sky_mean_scale(
            (self.center.x, self.center.y), wcs)
        inner_radius = self.inner_radius * mean_scale * u.arcsec
        outer_radius = self.outer_radius * mean_scale * u.arcsec
        return CircleAnnulusSkyRegion(center, inner_radius, outer_radius,
                                      self.meta.copy(), self.visual.copy())

    def to_polygon(self, *, n_vertices=100):
        """
        Return a `~regions.CompoundPixelRegion` of two
        `~regions.PolygonPixelRegion` objects that approximates this
        annulus.

        Parameters
        ----------
        n_vertices : int, optional
            The number of polygon vertices for each circle. Default
            is 100.

        Returns
        -------
        polygon : `~regions.CompoundPixelRegion`
            A compound region of two polygon regions approximating the
            annulus.
        """
        inner_polygon = self._inner_region.to_polygon(n_vertices=n_vertices)
        outer_polygon = self._outer_region.to_polygon(n_vertices=n_vertices)
        return CompoundPixelRegion(inner_polygon, outer_polygon,
                                   operator.xor, self.meta.copy(),
                                   self.visual.copy())

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        self._validate_planar_spherical_transform(
            wcs, include_boundary_distortions)

        if include_boundary_distortions:
            # Requires planar to spherical projection (using WCS)
            # and discretization
            # Will require implementing discretization in pixel space
            # to get correct handling of distortions.
            raise NotImplementedError

            # ### Potential solution:
            # # Leverage polygon class to_spherical_sky() functionality without
            # # distortions, as the distortions were already
            # # computed in creating
            # # that polygon approximation
            # return self.discretize_boundary(
            #     n_vertices=n_vertices).to_spherical_sky(
            #     wcs=wcs, include_boundary_distortions=False
            # )

        return self.to_sky(wcs).to_spherical_sky()


class CircleAnnulusSkyRegion(SkyRegion):
    """
    A circular annulus in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_radius : `~astropy.units.Quantity`
        The inner radius of the annulus in angular units.
    outer_radius : `~astropy.units.Quantity`
        The outer radius of the annulus in angular units.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_radius', 'outer_radius')
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_radius = PositiveScalarAngle('The inner radius as a |Quantity| '
                                       'angle.')
    outer_radius = PositiveScalarAngle('The outer radius as a |Quantity| '
                                       'angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_radius, outer_radius, meta=None,
                 visual=None):
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_radius >= outer_radius:
            raise ValueError('outer_radius must be greater than inner_radius')

    def to_pixel(self, wcs, *, as_ellipse=False):
        """
        Return a pixel region from this sky region.

        Parameters
        ----------
        wcs : WCS object
            A world coordinate system (WCS) transformation that
            supports the `astropy shared interface for WCS
            <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_
            (e.g., `astropy.wcs.WCS`).

        as_ellipse : bool, optional
            If `True`, return an `~regions.EllipseAnnulusPixelRegion`
            instead of a `~regions.CircleAnnulusPixelRegion`. An ellipse
            annulus is generally a better approximation when the WCS has
            distortions or different pixel scales along different axes.
            Default is `False`.

        Returns
        -------
        region : `~regions.CircleAnnulusPixelRegion` or `~regions.EllipseAnnulusPixelRegion`
            The pixel region. An ellipse annulus is returned if
            ``as_ellipse`` is `True`.
        """
        if as_ellipse:
            center, scale_major, scale_minor, angle = sky_to_pixel_svd_scales(
                self.center, wcs)
            inner_radius_arcsec = self.inner_radius.to_value(u.arcsec)
            outer_radius_arcsec = self.outer_radius.to_value(u.arcsec)
            inner_width = 2 * inner_radius_arcsec * scale_major
            outer_width = 2 * outer_radius_arcsec * scale_major
            inner_height = 2 * inner_radius_arcsec * scale_minor
            outer_height = 2 * outer_radius_arcsec * scale_minor
            return EllipseAnnulusPixelRegion(
                PixCoord(*center), inner_width, outer_width,
                inner_height, outer_height, angle=angle,
                meta=self.meta.copy(), visual=self.visual.copy())

        center, mean_scale = sky_to_pixel_mean_scale(self.center, wcs)
        inner_radius = self.inner_radius.to_value(u.arcsec) * mean_scale
        outer_radius = self.outer_radius.to_value(u.arcsec) * mean_scale
        return CircleAnnulusPixelRegion(PixCoord(*center), inner_radius,
                                        outer_radius,
                                        meta=self.meta.copy(),
                                        visual=self.visual.copy())

    def to_polygon(self, wcs, *, n_vertices=100):
        """
        Return a `~regions.CompoundSkyRegion` of two
        `~regions.PolygonSkyRegion` objects that approximates this
        annulus.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to use for the sky-to-pixel-to-sky conversion.
        n_vertices : int, optional
            The number of polygon vertices for each circle. Default
            is 100.

        Returns
        -------
        polygon : `~regions.CompoundSkyRegion`
            A compound region of two polygon regions approximating the
            annulus.
        """
        return self.to_pixel(wcs).to_polygon(n_vertices=n_vertices).to_sky(wcs)

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        self._validate_planar_spherical_transform(
            wcs, include_boundary_distortions)

        if include_boundary_distortions:
            # Requires planar to spherical projection (using WCS)
            # and discretization
            # Will require implementing discretization in pixel space
            # to get correct handling of distortions.
            raise NotImplementedError

            # ### Potential solution:
            # # Leverage polygon class to_spherical_sky() functionality without
            # # distortions, as the distortions were already
            # # computed in creating
            # # that polygon approximation
            # return self.to_pixel(wcs)
            #     .discretize_boundary(n_vertices=n_vertices)
            #     .to_spherical_sky(
            #     wcs=wcs, include_boundary_distortions=False
            # )

        return CircleAnnulusSphericalSkyRegion(
            self.center, self.inner_radius, self.outer_radius,
            self.meta.copy(), self.visual.copy(),
        )


class CircleAnnulusSphericalSkyRegion(AnnulusSphericalSkyRegion):
    """
    Class for a circular annulus sky region, where the circular annulus
    is interpreted within a spherical geometry reference frame.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The center position.
    inner_radius : `~astropy.units.Quantity`
        The inner radius in angular units.
    outer_radius : `~astropy.units.Quantity`
        The outer radius in angular units.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _component_class = CircleSphericalSkyRegion
    _params = ('center', 'inner_radius', 'outer_radius')
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_radius = PositiveScalarAngle('The inner radius as a |Quantity| '
                                       'angle.')
    outer_radius = PositiveScalarAngle('The outer radius as a |Quantity| '
                                       'angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_radius, outer_radius, meta=None,
                 visual=None):
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_radius >= outer_radius:
            raise ValueError('outer_radius must be greater than inner_radius')

    @property
    def _inner_region(self):
        return self._component_class(self.center, self.inner_radius,
                                     self.meta, self.visual)

    @property
    def _outer_region(self):
        return self._component_class(self.center, self.outer_radius,
                                     self.meta, self.visual)

    def transform_to(self, frame, merge_attributes=True):
        frame = self._validate_frame_transformation(frame)

        # Only center transforms, radii preserved
        center_transf = self.center.transform_to(
            frame, merge_attributes=merge_attributes)

        return CircleAnnulusSphericalSkyRegion(
            center_transf,
            self.inner_radius.copy(),
            self.outer_radius.copy(),
            self.meta.copy(),
            self.visual.copy(),
        )

    def discretize_boundary(self, *, n_vertices=100):
        return CompoundSphericalSkyRegion(
            self._inner_region.discretize_boundary(n_vertices=n_vertices),
            self._outer_region.discretize_boundary(n_vertices=n_vertices),
            operator=operator.xor,
            meta=self.meta.copy(),
            visual=self.visual.copy(),
        )

    def to_polygon(self, *, n_vertices=100):
        """
        Return a `~regions.CompoundSphericalSkyRegion` of two
        `~regions.PolygonSphericalSkyRegion` objects that approximates this
        annulus.

        Parameters
        ----------
        n_vertices : int, optional
            The number of polygon vertices for each circle. Default
            is 100.

        Returns
        -------
        polygon : `~regions.CompoundSphericalSkyRegion`
            A compound region of two polygon regions approximating the
            annulus.
        """
        return self.discretize_boundary(n_vertices=n_vertices)

    def to_sky(self, wcs=None, *, include_boundary_distortions=False,
               n_vertices=None):
        self._validate_planar_spherical_transform(
            wcs, include_boundary_distortions)

        if include_boundary_distortions:
            # Requires spherical to planar projection (from WCS)
            # and discretization
            # Use to_pixel(), then apply "small angle approx" to get planar
            # sky.
            return self.to_pixel(
                include_boundary_distortions=include_boundary_distortions,
                wcs=wcs, n_vertices=n_vertices,
            ).to_sky(wcs)

        return CircleAnnulusSkyRegion(
            self.center.copy(),
            self.inner_radius.copy(),
            self.outer_radius.copy(),
            meta=self.meta.copy(),
            visual=self.visual.copy(),
        )

    def to_pixel(self, wcs, *, include_boundary_distortions=False,
                 n_vertices=None):
        self._validate_planar_spherical_transform(
            wcs, include_boundary_distortions)

        if include_boundary_distortions:
            from .polygon import PolygonPixelRegion

            # Requires spherical to planar projection (from WCS) and
            # discretization
            disc_kwargs = (
                {} if n_vertices is None
                else {'n_vertices': n_vertices})
            polygonized = self.discretize_boundary(**disc_kwargs)

            inner_vertices = wcs.world_to_pixel(polygonized.region1.vertices)
            outer_vertices = wcs.world_to_pixel(polygonized.region2.vertices)

            return CompoundPixelRegion(
                PolygonPixelRegion(PixCoord(*inner_vertices)),
                PolygonPixelRegion(PixCoord(*outer_vertices)),
                operator=operator.xor,
                meta=self.meta.copy(),
                visual=self.visual.copy(),
            )

        return self.to_sky().to_pixel(wcs)


class AsymmetricAnnulusPixelRegion(AnnulusPixelRegion):
    """
    Helper class for asymmetric annuli sky regions.

    Used for ellipse and rectangle pixel annulus regions.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    inner_width : float
        The inner width of the annulus (before rotation) in pixels.
    outer_width : float
        The outer width of the annulus (before rotation) in pixels.
    inner_height : float
        The inner height of the annulus (before rotation) in pixels.
    outer_height : float
        The outer height of the annulus (before rotation) in pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the annulus, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the x
        axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_width', 'outer_width', 'inner_height',
               'outer_height', 'angle')
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        self.center = center
        self.inner_width = inner_width
        self.outer_width = outer_width
        self.inner_height = inner_height
        self.outer_height = outer_height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_width >= outer_width:
            raise ValueError('outer_width must be greater than inner_width')
        if inner_height >= outer_height:
            raise ValueError('outer_height must be greater than inner_height')

    @property
    def _inner_region(self):
        return self._component_class(self.center, self.inner_width,
                                     self.inner_height, self.angle,
                                     self.meta, self.visual)

    @property
    def _outer_region(self):
        return self._component_class(self.center, self.outer_width,
                                     self.outer_height, self.angle,
                                     self.meta, self.visual)

    def to_sky_args(self, wcs):
        # The photutils helpers measure the sky rotation as a position
        # angle (PA) from North; regions measures it from the RA axis.
        # Convert between them with a 90 deg offset.
        center, outer_width, outer_height, angle = pixel_shape_to_sky_svd(
            (self.center.x, self.center.y), wcs, self.outer_width,
            self.outer_height, self.angle.to_value(u.radian))
        _, inner_width, inner_height, _ = pixel_shape_to_sky_svd(
            (self.center.x, self.center.y), wcs, self.inner_width,
            self.inner_height, self.angle.to_value(u.radian))
        angle = (angle + 90 * u.deg).wrap_at(360 * u.deg)
        return (center, inner_width * u.arcsec, outer_width * u.arcsec,
                inner_height * u.arcsec, outer_height * u.arcsec, angle)


class AsymmetricAnnulusSkyRegion(SkyRegion):
    """
    Helper class for asymmetric annuli sky regions.

    Used for ellipse and rectangle sky annulus regions.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the annulus (before rotation) as an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the annulus (before rotation) as an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the annulus (before rotation) as an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the annulus (before rotation) as an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the annulus, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the
        longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_width', 'outer_width', 'inner_height',
               'outer_height', 'angle')

    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a'
                        '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        self.center = center
        self.inner_width = inner_width
        self.outer_width = outer_width
        self.inner_height = inner_height
        self.outer_height = outer_height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_width >= outer_width:
            raise ValueError('outer_width must be greater than inner_width')
        if inner_height >= outer_height:
            raise ValueError('outer_height must be greater than inner_height')

    def to_pixel_args(self, wcs):
        # Convert regions sky angle (from RA axis) to photutils PA (from
        # North) by subtracting 90 deg.
        sky_angle_rad = self.angle.to_value(u.radian) - math.pi / 2
        center, outer_width, outer_height, angle = sky_shape_to_pixel_svd(
            self.center, wcs,
            self.outer_width.to_value(u.arcsec),
            self.outer_height.to_value(u.arcsec),
            sky_angle_rad)
        _, inner_width, inner_height, _ = sky_shape_to_pixel_svd(
            self.center, wcs,
            self.inner_width.to_value(u.arcsec),
            self.inner_height.to_value(u.arcsec),
            sky_angle_rad)
        return (PixCoord(*center), inner_width, outer_width, inner_height,
                outer_height, angle)


class EllipseAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A elliptical annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the elliptical annulus.
    inner_width : float
        The inner width of the elliptical annulus (before rotation) in
        pixels.
    outer_width : float
        The outer width of the elliptical annulus (before rotation) in
        pixels.
    inner_height : float
        The inner height of the elliptical annulus (before rotation) in
        pixels.
    outer_height : float
        The outer height of the elliptical annulus (before rotation) in
        pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the elliptical annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the x axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.coordinates import Angle
        from regions import EllipseAnnulusPixelRegion, PixCoord

        fig, ax = plt.subplots()

        reg = EllipseAnnulusPixelRegion(PixCoord(6, 6),
                                        inner_width=5.5,
                                        outer_width=8.5,
                                        inner_height=3.5,
                                        outer_height=6.5,
                                        angle=Angle('45deg'))
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Ellipse Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    # duplicated from AsymmetricAnnulusPixelRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = EllipsePixelRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_polygon(self, *, n_vertices=100):
        """
        Return a `~regions.CompoundPixelRegion` of two
        `~regions.PolygonPixelRegion` objects that approximates this
        annulus.

        Parameters
        ----------
        n_vertices : int, optional
            The number of polygon vertices for each ellipse. Default
            is 100.

        Returns
        -------
        polygon : `~regions.CompoundPixelRegion`
            A compound region of two polygon regions approximating the
            annulus.
        """
        inner_polygon = self._inner_region.to_polygon(n_vertices=n_vertices)
        outer_polygon = self._outer_region.to_polygon(n_vertices=n_vertices)
        return CompoundPixelRegion(inner_polygon, outer_polygon,
                                   operator.xor, self.meta.copy(),
                                   self.visual.copy())

    def to_sky(self, wcs):
        return EllipseAnnulusSkyRegion(*self.to_sky_args(wcs),
                                       meta=self.meta.copy(),
                                       visual=self.visual.copy())

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        raise NotImplementedError


class EllipseAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A elliptical annulus in `~astropy.coordinates.SkyCoord` coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the elliptical annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the elliptical annulus (before rotation) as
        an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the elliptical annulus (before rotation) as
        an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the elliptical annulus (before rotation) as
        an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the elliptical annulus (before rotation) as
        an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the elliptical annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    # duplicated from AsymmetricAnnulusSkyRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = EllipseSkyRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_polygon(self, wcs, *, n_vertices=100):
        """
        Return a `~regions.CompoundSkyRegion` of two
        `~regions.PolygonSkyRegion` objects that approximates this
        annulus.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to use for the sky-to-pixel-to-sky conversion.
        n_vertices : int, optional
            The number of polygon vertices for each ellipse. Default
            is 100.

        Returns
        -------
        polygon : `~regions.CompoundSkyRegion`
            A compound region of two polygon regions approximating the
            annulus.
        """
        return self.to_pixel(wcs).to_polygon(n_vertices=n_vertices).to_sky(wcs)

    def to_pixel(self, wcs):
        return EllipseAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                         meta=self.meta.copy(),
                                         visual=self.visual.copy())

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        raise NotImplementedError


class RectangleAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A rectangular annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the rectangular annulus.
    inner_width : float
        The inner width of the rectangular annulus (before rotation) in
        pixels.
    outer_width : float
        The outer width of the rectangular annulus (before rotation) in
        pixels.
    inner_height : float
        The inner height of the rectangular annulus (before rotation) in
        pixels.
    outer_height : float
        The outer height of the rectangular annulus (before rotation) in
        pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangular annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the x axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from astropy.coordinates import Angle
        from regions import PixCoord, RectangleAnnulusPixelRegion

        fig, ax = plt.subplots()

        reg = RectangleAnnulusPixelRegion(PixCoord(x=6, y=6),
                                          inner_width=5.5,
                                          outer_width=8.5,
                                          inner_height=3.5,
                                          outer_height=6.5,
                                          angle=Angle('45deg'))
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Rectangle Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    # duplicated from AsymmetricAnnulusPixelRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = RectanglePixelRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_sky_args(self, wcs):
        # The photutils SVD helpers measure the sky rotation as a
        # position angle (PA) from North; regions measures it from the
        # RA axis. Convert between them with a 90 deg offset.
        center, outer_width, outer_height, angle = pixel_shape_to_sky_svd(
            (self.center.x, self.center.y), wcs, self.outer_width,
            self.outer_height, self.angle.to_value(u.radian))
        _, inner_width, inner_height, _ = pixel_shape_to_sky_svd(
            (self.center.x, self.center.y), wcs, self.inner_width,
            self.inner_height, self.angle.to_value(u.radian))
        angle = (angle + 90 * u.deg).wrap_at(360 * u.deg)
        return (center, inner_width * u.arcsec, outer_width * u.arcsec,
                inner_height * u.arcsec, outer_height * u.arcsec, angle)

    def to_polygon(self):
        """
        Return a `~regions.CompoundPixelRegion` of two
        `~regions.PolygonPixelRegion` objects equivalent to this
        annulus.

        Returns
        -------
        polygon : `~regions.CompoundPixelRegion`
            A compound region of two polygon regions equivalent to the
            annulus.
        """
        inner_polygon = self._inner_region.to_polygon()
        outer_polygon = self._outer_region.to_polygon()
        return CompoundPixelRegion(inner_polygon, outer_polygon,
                                   operator.xor, self.meta.copy(),
                                   self.visual.copy())

    def to_sky(self, wcs):
        return RectangleAnnulusSkyRegion(*self.to_sky_args(wcs),
                                         meta=self.meta.copy(),
                                         visual=self.visual.copy())

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        raise NotImplementedError


class RectangleAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A rectangular annulus in `~astropy.coordinates.SkyCoord`
    coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangular annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the rectangular annulus (before rotation) as
        an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the rectangular annulus (before rotation) as
        an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the rectangular annulus (before rotation) as
        an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the rectangular annulus (before rotation) as
        an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangular annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    # duplicated from AsymmetricAnnulusSkyRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = RectangleSkyRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_pixel_args(self, wcs):
        # Convert regions sky angle (from RA axis) to photutils PA (from
        # North) by subtracting 90 deg.
        sky_angle_rad = self.angle.to_value(u.radian) - math.pi / 2
        center, outer_width, outer_height, angle = sky_shape_to_pixel_svd(
            self.center, wcs,
            self.outer_width.to_value(u.arcsec),
            self.outer_height.to_value(u.arcsec),
            sky_angle_rad)
        _, inner_width, inner_height, _ = sky_shape_to_pixel_svd(
            self.center, wcs,
            self.inner_width.to_value(u.arcsec),
            self.inner_height.to_value(u.arcsec),
            sky_angle_rad)
        return (PixCoord(*center), inner_width, outer_width, inner_height,
                outer_height, angle)

    def to_polygon(self, wcs):
        """
        Return a `~regions.CompoundSkyRegion` of two
        `~regions.PolygonSkyRegion` objects equivalent to this
        annulus.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to use for the sky-to-pixel-to-sky conversion.

        Returns
        -------
        polygon : `~regions.CompoundSkyRegion`
            A compound region of two polygon regions equivalent to the
            annulus.
        """
        return self.to_pixel(wcs).to_polygon().to_sky(wcs)

    def to_pixel(self, wcs):
        return RectangleAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                           meta=self.meta.copy(),
                                           visual=self.visual.copy())

    def to_spherical_sky(self, wcs=None, *, include_boundary_distortions=False,
                         n_vertices=None):
        raise NotImplementedError
