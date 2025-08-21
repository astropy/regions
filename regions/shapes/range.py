# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines a range spherical sky region (bounded by longitude
and/or latitude ranges).
"""

import copy
import operator

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.stats import circmean

from regions._utils.spherical_helpers import (
    cross_product_skycoord2skycoord, cross_product_sum_skycoord2skycoord,
    discretize_all_edge_boundaries, get_edge_raw_lonlat_bounds_circ_edges)
from regions.core.attributes import TwoValAngleorNone
from regions.core.compound import CompoundSphericalSkyRegion
from regions.core.core import ComplexSphericalSkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord
from regions.shapes.annulus import CircleAnnulusSphericalSkyRegion
from regions.shapes.circle import CircleSphericalSkyRegion
from regions.shapes.lune import LuneSphericalSkyRegion
from regions.shapes.polygon import (PolygonPixelRegion, PolygonSkyRegion,
                                    PolygonSphericalSkyRegion)
from regions.shapes.whole_sky import WholeSphericalSkyRegion

__all__ = ['RangeSphericalSkyRegion']


class RangeSphericalSkyRegion(ComplexSphericalSkyRegion):
    """
    Sky Region defined within a range of longitude and/or latitudes. At
    least one set of longitude or latitude bounds must be set.

    If latitude_range[0] > latitude_range[1], will wrap over the poles
    (eg, will *exclude* the central latitude range.)

    Parameters
    ----------
    frame : `~astropy.coordinates.BaseCoordinateFrame`
        The coordinate frame for the specified range region.

    longitude_range : length 2 list-like of `~astropy.coordinates.Angle` or
        `~astropy.unit.Quantity` or None
        Longitude range in region. Spans from first to second entries
        (in increasing longitude value), so inverting the order will select
        the complement longitude range.

    latitude_range : length 2 list-like of `~astropy.coordinates.Angle` or
        `~astropy.unit.Quantity` or None
        Latitude range in region. Spans from first to second entries
        (in increasing latitude value), so inverting the order will select
        the complement latitude range.
        I.e., if latitude_range[0] > latitude_range[1], will wrap over the poles,
        and will instead exclude the central latitude range.

    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.

    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    **kwargs : Keyword arguments
        Additional keyword arguments passed when created new instances
        after coordinate transformations.
    """

    _params = ('frame', 'longitude_range', 'latitude_range')
    _boundaries = ('longitude_bounds', 'latitude_bounds')
    longitude_range = TwoValAngleorNone(
        'The range of longitude values, as a length-2 |Quantity| angle or list or as None'
    )
    latitude_range = TwoValAngleorNone(
        'The range of longitude values, as a length-2 |Quantity| angle or list or as None'
    )

    def __init__(
        self,
        frame=None,
        longitude_range=None,
        latitude_range=None,
        meta=None,
        visual=None,
        **kwargs
    ):

        self.longitude_range = longitude_range
        self.latitude_range = latitude_range

        self._longitude_bounds = None
        self._latitude_bounds = None
        self._vertices = None
        self._is_original_frame = kwargs.pop('_is_original_frame', True)

        # TODO:
        # Validate lon/lat range inputs, to ensure lat range has expected definition.

        # Check if "_frame" is passed as a kwarg, from a copy:
        _frame = None
        if kwargs.get('_frame', 'missing') != 'missing':
            _frame = kwargs.pop('_frame', None)

        # Passing frame overrides any value passed in "_frame" in the kwargs,
        # so only if "frame"=None is "_frame" used:
        if frame is None:
            frame = _frame

        # Validate frame and get into standard format:
        frame = self._validate_frame(frame)

        # --------------------
        if (
            (longitude_range is None)
            & (latitude_range is None)
        ):
            if (
                (kwargs.get('_longitude_bounds', None) is not None)
                | (kwargs.get('_latitude_bounds', None) is not None)
            ):
                # Create class by directly passing long/lat bounds
                # -- for cases when creating transformed RangeSphericalSkyRegion instances

                # TODO: validation of direct bounds inputs
                self._longitude_bounds = kwargs.pop('_longitude_bounds', None)
                self._latitude_bounds = kwargs.pop('_latitude_bounds', None)

                # Also directly input _vertices, as these are only computed in the
                # original frame
                self._vertices = kwargs.pop('_vertices', None)

                # Also set _params to to just ("frame") for this transformed instance.
                self._params = ('frame',)

            else:
                raise ValueError(
                    'A range for at least one of longitude and latitude must be set'
                )

        # Stash frame in private attribute to access when building on-the-fly
        # derived attributes & bounds:
        self._frame = frame

        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def __eq__(self, other):
        # Equality operator for RangeSphericalSkyRegion.

        # Modified from default, because the lon/lat ranges can be none
        # -- requiring comparison of boundaries.

        # All Region properties are compared for strict equality except
        # for Quantity parameters, which allow for different units if they
        # are directly convertible.

        if not isinstance(other, self.__class__):
            return False

        meta_params = ['meta', 'visual']
        intern_params = ['_frame', '_vertices', '_is_original_frame', '_params']
        if self._params is not None:
            self_params = (
                list(self._params) + [f"_{bn}" for bn in list(self._boundaries)]
                + intern_params + meta_params
            )
        else:
            self_params = (
                [f"_{bn}" for bn in list(self._boundaries)]
                + intern_params + meta_params
            )

        if other._params is not None:
            other_params = (
                list(other._params) + [f"_{bn}" for bn in list(other._boundaries)]
                + intern_params + meta_params
            )
        else:
            other_params = (
                [f"_{bn}" for bn in list(other._boundaries)]
                + intern_params + meta_params
            )

        # check that both have identical parameters
        if self_params != other_params:
            return False

        # now check the parameter values
        # Note that Quantity comparisons allow for different units
        # if they directly convertible (e.g., 1. * u.deg == 60. * u.arcmin)
        try:
            for param in self_params:
                # np.any is used for SkyCoord array comparisons
                if np.any(getattr(self, param) != getattr(other, param)):
                    return False
        except TypeError:
            # TypeError is raised from SkyCoord comparison when they do
            # not have equivalent frames. Here return False instead of
            # the TypeError.
            return False

        return True

    # ALWAYS derive boundaries on the fly, IF all _params are not None
    # --- only a concern for RangeSphericalSkyRegion....

    @property
    def __longitude_bounds(self):
        # Internal, on-the-fly boundaries determination.
        # A concern if range values change after original _longitude_bounds
        # are computed, if these were "static" for all cases
        if self.longitude_range is None:
            return None

        # Define GC centers based on longitude range
        c1 = SkyCoord(self.longitude_range[0], 0 * u.deg, frame=self.frame)
        c2 = SkyCoord(self.longitude_range[1], 0 * u.deg, frame=self.frame)

        # Cross product to get pole.
        # Invert order if lon_range[1] > lon_range[0]
        if self.longitude_range[0] > self.longitude_range[1]:
            c_pole = cross_product_skycoord2skycoord(c2, c1)
        else:
            c_pole = cross_product_skycoord2skycoord(c1, c2)

        # GC centers from cross products:
        center_gc1 = cross_product_skycoord2skycoord(c_pole, c1)
        center_gc2 = cross_product_skycoord2skycoord(c2, c_pole)

        c_gc1_sph = center_gc1.represent_as('spherical')
        c_gc2_sph = center_gc2.represent_as('spherical')

        # Squash small machine-rounding problems in cross product:
        # latitude should be 0 for these longitude-describing GCs,
        # as they pass through the poles.
        # Only possible in the original frame, when defining from longitude_range

        center_gc1 = SkyCoord(c_gc1_sph.lon, 0 * u.deg, frame=self.frame)
        center_gc2 = SkyCoord(c_gc2_sph.lon, 0 * u.deg, frame=self.frame)

        return LuneSphericalSkyRegion(
            center_gc1, center_gc2, self.meta, self.visual
        )

    @property
    def __latitude_bounds(self):
        # Internal, on-the-fly boundaries determination.
        # A concern if range values change after original _latitude_bounds
        # are computed, if these were "static" for all cases
        if self.latitude_range is None:
            return None

        # Define north pole: 0 lon, 90 lat
        c_pole = SkyCoord(0 * u.deg, 90 * u.deg, frame=self.frame)

        lat_arr = self.latitude_range.copy()

        # Logic if this is a "across-polar" region is handled by swapping
        # the logic of the latitude boundary itself:
        _wraps_poles = False
        if lat_arr[0] > lat_arr[1]:
            _wraps_poles = True
            lat_arr = lat_arr[::-1]

        inner_rad = 90.0 * u.deg - lat_arr[1]
        outer_rad = 90.0 * u.deg - lat_arr[0]

        # Special handling: check if either edge of the range is -90 or 90:
        if (
            (lat_arr[1] == Angle(90 * u.deg))
            & (lat_arr[0] == Angle(-90 * u.deg))
        ):
            # Whole sky: no constraint:
            return None

        elif lat_arr[1] == Angle(90 * u.deg):
            if _wraps_poles:
                # Check if it wrapped pole, with input [90*u.deg, X*u.deg]:
                # equivalent to [-90*u.deg, X*u.deg], and is a negated simple circle
                # Radius from input X*u.deg constraint is in outer_rad after swap
                CompoundSphericalSkyRegion(
                    WholeSphericalSkyRegion(),
                    CircleSphericalSkyRegion(
                        c_pole, outer_rad, self.meta, self.visual
                    ),
                    operator.xor, self.meta, self.visual
                )
            else:
                # Simple circle case:
                return CircleSphericalSkyRegion(
                    c_pole, outer_rad, self.meta, self.visual
                )
        elif lat_arr[0] == Angle(-90 * u.deg):
            if _wraps_poles:
                # Check if it wrapped pole, with input [X*u.deg, -90*u.deg],
                # equivalent to [X*u.deg, 90*u.deg], and is a simple circle
                # Radius from input X*u.deg constraint is in inner_rad after swap
                CircleSphericalSkyRegion(
                    c_pole, inner_rad, self.meta, self.visual
                )
            else:
                # Negated simple circle
                return CompoundSphericalSkyRegion(
                    WholeSphericalSkyRegion(),
                    CircleSphericalSkyRegion(
                        c_pole, inner_rad, self.meta, self.visual
                    ),
                    operator.xor, self.meta, self.visual
                )

        # True annular range:
        if _wraps_poles:
            # Negate by passing as xor with whole sky
            return CompoundSphericalSkyRegion(
                WholeSphericalSkyRegion(),
                CircleAnnulusSphericalSkyRegion(
                    c_pole, inner_rad, outer_rad,
                    self.meta, self.visual
                ),
                operator.xor, self.meta, self.visual
            )

        # Not wrapping poles:
        return CircleAnnulusSphericalSkyRegion(
            c_pole, inner_rad, outer_rad, self.meta, self.visual
        )

    @property
    def _bound_nverts(self):
        # Number of vertices for boundary:
        # lat-only: 0 (annulus)
        # lon-only: 2 (lune)
        # both: 4 (range)

        # Only compute params once:
        lon_bounds = self.longitude_bounds
        lat_bounds = self.latitude_bounds

        if (lon_bounds is not None) & (lat_bounds is not None):
            # 4 vertices:
            nverts = 4
        elif (lon_bounds is not None) & (lat_bounds is None):
            # 2 vertices:
            nverts = 2
        else:
            # Only lat_bounds: 0 vertices
            nverts = 0
        return nverts

    @property
    def _edge_circs(self):
        """
        Get list of the circles defining the range boundaries.

        Only used if vertices=4, otherwise lune/annulus properties used.
        """
        bound_nverts = self._bound_nverts

        if bound_nverts == 4:
            return [
                self.longitude_bounds._circle_1,
                self.latitude_bounds._outer_region,
                self.longitude_bounds._circle_2,
                self.latitude_bounds._inner_region,
            ]

        # Otherwise, return None
        return None

    @property
    def _compound_region(self):
        if self.longitude_bounds is None:
            return self.latitude_bounds

        if self.latitude_bounds is None:
            return self.longitude_bounds

        return CompoundSphericalSkyRegion(
            self.longitude_bounds, self.latitude_bounds,
            operator.and_, self.meta, self.visual
        )

    @property
    def frame(self):
        return self._frame

    @property
    def longitude_bounds(self):
        """
        Longitude bounds of the range spherical sky region.

        Defined as a spherical lune, or None if no longitude bound set.
        """
        # If only boundaries are set: get bounds from internal attribute:
        if (
            (self.longitude_range is not None)
            & (getattr(self, '_longitude_bounds', None) is not None)
        ):
            raise ValueError('Both longitude_range and _longitude_bounds set!')

        if getattr(self, '_longitude_bounds', None) is not None:
            return self._longitude_bounds

        # Otherwise derive on-the-fly
        return self.__longitude_bounds

    @property
    def latitude_bounds(self):
        """
        Latitude bounds of the range spherical sky region.

        Defined as a spherical circular annulus, a circle (if ending at
        a pole), or None if no latitude bound set.
        """
        # If only boundaries are set: get bounds from internal attribute:
        if (
            (self.latitude_range is not None)
            & (getattr(self, '_latitude_bounds', None) is not None)
        ):
            raise ValueError('Both latitude_range and _latitude_bounds set!')

        if getattr(self, '_latitude_bounds', None) is not None:
            return self._latitude_bounds

        # Otherwise derive on-the-fly
        return self.__latitude_bounds

    @property
    def centroid(self):
        """
        Region centroid.

        If both longitude and latitude bounds are set:

        Defined as the point equidistant from all vertices.

        However, if this point falls outside of the region
        (as in cases very "long and narrow" ranges),
        the centroid is instead approximated as the average
        of the longitude and latitudes of all vertices.

        If only longitude or only latitude bounds are set,
        the centroid is taken as the centroid of that
        boundary shape (i.e., the lune centroid or annulus center).
        """
        return self.centroid_mindist

    @property
    def centroid_avg(self):
        """
        An approximate "centroid", taking the average of the vertices
        longitude / latitudes, and invert as needed to define an
        "inside" point for a wide angle longitude range.

        Will behave poorly for regions crossing the poles.
        """
        bound_nverts = self._bound_nverts

        if bound_nverts == 4:
            # Both lon and lat bounds
            verts_sph = self.vertices.represent_as('spherical')

            lons = verts_sph.lon
            lats = verts_sph.lat

            lon = circmean(lons)
            lat = np.average(lats)
            centroid = SkyCoord(lon, lat, frame=self.frame)
            if (
                (not self.contains(centroid)) and (
                    not self.longitude_bounds.contains(centroid)
                )
            ):
                # If centroid not contained in range, or in subset longitude bounds:
                # modify longitude by adding 180 deg:
                centroid = SkyCoord(lon + 180 * u.deg, lat, frame=self.frame)

            return centroid
        elif bound_nverts == 2:
            # Only lon bounds: just get centroid of lune:
            return self.longitude_bounds.centroid
        elif bound_nverts == 0:
            # Only lat_bounds: just get center of annulus:
            return self.latitude_bounds.center

    @property
    def centroid_mindist(self):
        """
        Region centroid.

        Defined as the point equidistant from all vertices.
        """
        bound_nverts = self._bound_nverts

        if bound_nverts == 4:
            # If both lon/lat bounds:
            # Calculate from sum of cross products of vertices with cartesian representation
            # verts are in CW order:
            centroid_mindist = cross_product_sum_skycoord2skycoord(self.vertices)

            if not self.contains(centroid_mindist):
                # If it's not contained, check if it's because it's
                # outside the latitude range: if so, instead just use avg centroid:
                if not self.latitude_bounds.contains(centroid_mindist):
                    # Problem of very elongated bounds at high latitude --
                    # end up outside range. If so, fall back to average
                    return self.centroid_avg
                raise ValueError

            return centroid_mindist
        elif bound_nverts == 2:
            # Only lon bounds: just get centroid of lune:
            return self.longitude_bounds.centroid
        elif bound_nverts == 0:
            # Only lat_bounds: just get center of annulus:
            return self.latitude_bounds.center

    @property
    def vertices(self):
        """
        Region vertices.

        The result depends on which ranges are set:

        If both longitude and latitude ranges are set, returns 4 vertices
        (the "corners" of the composite range) as a SkyCoord.

        If only a longitude range is set, returns 2 vertices (the
        vertices of the longitude spherical lune boundary) as a SkyCoord.

        If only latitude ranges are set, returns None (no vertices
        for the spherical circular annulus latitude bounds).
        """
        if getattr(self, '_vertices', None) is not None:
            return self._vertices
        else:
            bound_nverts = self._bound_nverts

            if bound_nverts == 4:
                # Both lon/lat bounds:on the fly construct from the lon/lat range attrib
                # Longitude min is on the right, so derive CW from lower right corner.

                # Note this is never called for transformed instances,
                # which have vertices info stashed in _vertices.
                verts_lon = [
                    self.longitude_range[0],
                    self.longitude_range[1],
                    self.longitude_range[1],
                    self.longitude_range[0]
                ]
                verts_lat = [
                    self.latitude_range[0],
                    self.latitude_range[0],
                    self.latitude_range[1],
                    self.latitude_range[1],
                ]
                return SkyCoord(verts_lon, verts_lat, frame=self.frame)

            if bound_nverts == 2:
                # 2 vertices: directly return from lune itself.
                return self.longitude_bounds.vertices
            elif bound_nverts == 0:
                # Only lat_bounds: no vertices
                return None

    @property
    def bounding_circle(self):
        bound_nverts = self._bound_nverts
        if bound_nverts == 0:
            # No vertices: only an annulus:
            # Use the annulus bounding circle:
            return self.latitude_bounds.bounding_circle()
        elif bound_nverts == 2:
            # Lune: get lune bounding circle
            return self.longitude_bounds.bounding_circle()
        elif bound_nverts == 4:
            # Quadrilateral: use max of seps to vertices:
            cent = copy.deepcopy(self.centroid_mindist)
            seps = cent.separation(self.vertices)
            return CircleSphericalSkyRegion(center=cent, radius=np.max(seps))
        else:
            raise ValueError

    @property
    def bounding_lonlat(self):
        bound_nverts = self._bound_nverts
        if bound_nverts == 0:
            # Only an annulus: Use the annulus bounding lonlat:
            return self.latitude_bounds.bounding_lonlat()

        if bound_nverts == 2:
            # Lune: get lune bounding lonlat
            return self.longitude_bounds.bounding_lonlat()

        if bound_nverts == 4:
            # Only used for both lon+lat bounds,
            # otherwise directly calls lon or lat bounding_lonlat()
            lons_arr, lats_arr = get_edge_raw_lonlat_bounds_circ_edges(
                self.vertices, self.centroid, self._edge_circs,
            )

            lons_arr, lats_arr = self._validate_lonlat_bounds(lons_arr, lats_arr)

            return lons_arr, lats_arr

    def contains(self, coord):
        if self._is_original_frame:
            # Just compute directly from lon/lat ranges:

            # Ensure coordinate is in same frame:
            c_repr = coord.represent_as('spherical')

            # Longitude range constraints:
            if self.longitude_range is not None:
                # Ensure lon range is wrapped at 360, matching internal SkyCoord convention
                # To handle possible cases of lon ranges "wrapping around" across
                # the standard 360->0 wrap, wrap lon ranges + coord values
                # around the first entry angle
                wrap_ang = self.longitude_range[0]
                coord_lons = c_repr.lon.wrap_at(wrap_ang)
                in_range_lon = (
                    (coord_lons >= Angle(self.longitude_range[0]).wrap_at(wrap_ang))
                    & (coord_lons <= Angle(self.longitude_range[1]).wrap_at(wrap_ang))
                )
            # No constraints
            elif coord.isscalar:
                in_range_lon = True
            else:
                in_range_lon = np.ones(coord.shape, dtype=bool)

            # Latitude range constraints:
            if self.latitude_range is not None:
                coord_lat = c_repr.lat

                # Check if latitude ranges are "simpler":
                # Equivalent to no constraints, or simpler circles:
                if (
                    (self.latitude_range[1] == Angle(90 * u.deg))
                    & (self.latitude_range[0] == Angle(-90 * u.deg))
                ):
                    # Whole sky:
                    if coord.isscalar:
                        in_range_lat = True
                    else:
                        in_range_lat = np.ones(coord.shape, dtype=bool)
                elif (
                    self.latitude_range[1] == Angle(90 * u.deg)
                ):
                    # Range includes latitude up to N pole
                    in_range_lat = coord_lat >= self.latitude_range[0]
                elif (
                    self.latitude_range[0] == Angle(-90 * u.deg)
                ):
                    # Range includes latitudes down to S pole
                    in_range_lat = coord_lat <= self.latitude_range[1]
                else:
                    # Actual annulus
                    in_range_lat = (
                        (coord_lat >= self.latitude_range[0])
                        & (coord_lat <= self.latitude_range[1])
                    )
            # No constraints
            elif coord.isscalar:
                in_range_lat = True
            else:
                in_range_lat = np.ones(coord.shape, dtype=bool)

            in_range = in_range_lon & in_range_lat

            if self.meta.get('include', True):
                return in_range
            else:
                return np.logical_not(in_range)

        # If transformed: must use boundary compound region logic:
        return self._compound_region.contains(coord)

    def transform_to(self, frame, merge_attributes=True):
        frame = self._validate_frame(frame)

        transformed = {}

        # Internal attribute noting this is a transformed frame.
        transformed['_is_original_frame'] = False

        fields = list(self._params)
        for field in fields:
            # High-level lon/lat attributes will not transform: set to None
            transformed[field] = None

        boundaries = list(self._boundaries)
        for bound in boundaries:
            # Will be passing directly into internal attributes
            # _longitude_bounds, _latitude_bounds
            transformed[f"_{bound}"] = getattr(self, bound).transform_to(
                frame, merge_attributes=merge_attributes
            )

        # Also transform vertices:
        for field in ['_vertices']:
            verts = copy.deepcopy(getattr(self, field))
            if (verts is None) and self._is_original_frame:
                # If internal attribute not set, and original frame,
                # get the on-the-fly-derived verts.
                verts = self.vertices
            if verts is not None:
                transformed[field] = verts.transform_to(
                    frame, merge_attributes=merge_attributes
                )
            else:
                transformed[field] = None

        # Add transformed frame to transformed dict:
        transformed['frame'] = frame

        return self.__class__(**transformed)

    def discretize_boundary(self, n_points=10):
        bound_nverts = self._bound_nverts

        if bound_nverts == 0:
            return self.latitude_bounds.discretize_boundary(n_points=n_points)
        elif bound_nverts == 2:
            return self.longitude_bounds.discretize_boundary(n_points=n_points)
        elif bound_nverts == 4:
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
            disc_bound = self.discretize_boundary(**discretize_kwargs)
            if isinstance(disc_bound, CompoundSphericalSkyRegion):
                return disc_bound.to_pixel(wcs)

            verts = wcs.world_to_pixel(disc_bound.vertices)

            return PolygonPixelRegion(
                PixCoord(*verts), meta=self.meta.copy(),
                visual=self.visual.copy()
            )

        return self.to_sky().to_pixel(wcs)
