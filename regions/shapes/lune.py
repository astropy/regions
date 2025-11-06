# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines a lune (a "wedge" bounded by two great circles)
spherical sky region.
"""

import operator

import astropy.units as u
from astropy.coordinates import SkyCoord

from regions._utils.spherical_helpers import (
    cross_product_skycoord2skycoord, cross_product_sum_skycoord2skycoord,
    discretize_all_edge_boundaries, get_edge_raw_lonlat_bounds_circ_edges)
from regions.core.attributes import (RegionMetaDescr, RegionVisualDescr,
                                     ScalarSkyCoord)
from regions.core.compound import CompoundSphericalSkyRegion
from regions.core.core import SphericalSkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.shapes.circle import CircleSphericalSkyRegion
from regions.shapes.polygon import PolygonSphericalSkyRegion

__all__ = ['LuneSphericalSkyRegion']


class LuneSphericalSkyRegion(SphericalSkyRegion):
    """
    Class for a spherical "lune", the intersection between to great
    circles in spherical geometry.

    Composed of two CircleSphericalSkyRegions.

    Parameters
    ----------
    center_gc1 : `~astropy.coordinates.SkyCoord`
        The center position of the first great circle.
    center_gc2 : `~astropy.coordinates.SkyCoord`
        The center position of the second great circle.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center_gc1', 'center_gc2')
    center_gc1 = ScalarSkyCoord('The center position of the first great circle as a |SkyCoord|.')
    center_gc2 = ScalarSkyCoord('The center position of the second great circle as a |SkyCoord|.')

    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center_gc1, center_gc2, meta=None, visual=None):
        self.center_gc1 = center_gc1
        self.center_gc2 = center_gc2

        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def _circle_1(self):
        return CircleSphericalSkyRegion(
            self.center_gc1, 90 * u.deg,
            self.meta, self.visual
        )

    @property
    def _circle_2(self):
        return CircleSphericalSkyRegion(
            self.center_gc2, 90 * u.deg,
            self.meta, self.visual
        )

    @property
    def _edge_circs(self):
        """
        Get list of the great circles defining the lune boundaries.
        """
        return [self._circle_1, self._circle_2]

    @property
    def _compound_region(self):
        return CompoundSphericalSkyRegion(self._circle_1, self._circle_2,
                                          operator.and_, self.meta, self.visual)

    @property
    def frame(self):
        return self.center_gc1.frame

    @property
    def centroid(self):
        """
        Region centroid.
        """
        # Cross product to get poles:
        c_pole = cross_product_skycoord2skycoord(self.center_gc2, self.center_gc1)
        c_pole_a = cross_product_skycoord2skycoord(self.center_gc1, self.center_gc2)

        # Use centers + poles to get centroid from cross products of vertices
        # with cartesian representation
        # Assume CW order is gc1, c_pole, gc2, c_pole_a

        c_gc1_sph = self.center_gc1.represent_as('spherical')
        c_gc2_sph = self.center_gc2.represent_as('spherical')
        c_p_sph = c_pole.represent_as('spherical')
        c_p_a_sph = c_pole_a.represent_as('spherical')

        verts = SkyCoord(
            [
                c_gc2_sph.lon,
                c_p_a_sph.lon,
                c_gc1_sph.lon,
                c_p_sph.lon,
            ],
            [
                c_gc2_sph.lat,
                c_p_a_sph.lat,
                c_gc1_sph.lat,
                c_p_sph.lat,
            ],
            frame=self.frame
        )

        return cross_product_sum_skycoord2skycoord(verts)

    @property
    def vertices(self):
        """
        Region vertices.

        For a lune, these are the intersections of the two bounding
        great circles, and are antipodes.
        """
        # Cross product to get poles:
        c_pole = cross_product_skycoord2skycoord(self.center_gc2, self.center_gc1)
        c_pole_a = cross_product_skycoord2skycoord(self.center_gc1, self.center_gc2)

        c_p_sph = c_pole.represent_as('spherical')
        c_p_a_sph = c_pole_a.represent_as('spherical')

        verts = SkyCoord(
            [
                c_p_sph.lon,
                c_p_a_sph.lon,
            ],
            [
                c_p_sph.lat,
                c_p_a_sph.lat,
            ],
            frame=self.frame
        )

        return verts

    @property
    def bounding_circle(self):
        return CircleSphericalSkyRegion(self.centroid, 90 * u.deg)

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

        center_gc1_transf = self.center_gc1.transform_to(frame, merge_attributes=merge_attributes)
        center_gc2_transf = self.center_gc2.transform_to(frame, merge_attributes=merge_attributes)

        return LuneSphericalSkyRegion(
            center_gc1_transf,
            center_gc2_transf,
            self.meta.copy(),
            self.visual.copy()
        )

    def discretize_boundary(self, n_points=100):
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
        if not include_boundary_distortions:
            raise ValueError(
                'Invalid parameter: `include_boundary_distortions=False`!\n'
                'Transforming lune to planar sky region is only possible by '
                'including boundary distortions, as there is no analogous sky region.'
            )
        if discretize_kwargs is None:
            discretize_kwargs = {}

        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires spherical to planar projection (from WCS) and discretization
            # Use to_pixel(), then apply "small angle approx" to get planar sky.
            return self.to_pixel(
                include_boundary_distortions=include_boundary_distortions,
                wcs=wcs,
                discretize_kwargs=discretize_kwargs,
            ).to_sky(wcs)

    def to_pixel(
            self,
            wcs=None,
            include_boundary_distortions=False,
            discretize_kwargs=None,
    ):
        if not include_boundary_distortions:
            raise ValueError(
                'Invalid parameter: `include_boundary_distortions=False`!\n'
                'Transforming range to planar pixel region is only possible by '
                'including boundary distortions, as there is no analogous pixel region.'
            )

        if discretize_kwargs is None:
            discretize_kwargs = {}
        if include_boundary_distortions:
            if wcs is None:
                raise ValueError(
                    "'wcs' must be set if 'include_boundary_distortions'=True"
                )
            # Requires spherical to planar projection (from WCS) and discretization
            raise NotImplementedError
