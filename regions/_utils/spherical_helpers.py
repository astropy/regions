# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides spherical sky region calculation helper tools.
"""

import astropy.units as u
import numpy as np
from astropy.coordinates import (Latitude, Longitude, SkyCoord,
                                 SphericalRepresentation,
                                 UnitSphericalRepresentation,
                                 cartesian_to_spherical, concatenate)

__all__ = []


def cross_product_skycoord2skycoord(c1, c2):
    """
    Compute cross product of two sky coordinates (from a spherical
    representation), returning a third sky coordinate (on a spherical
    representation).

    Parameters
    ----------
    c1 : `~astropy.coordinates.SkyCoord`
        The first SkyCoord.

    c2 : `~astropy.coordinates.SkyCoord`
        The second SkyCoord.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The cross product as a SkyCoord, with a spherical representation.
    """
    cross = c1.frame.represent_as('cartesian').cross(c2.frame.represent_as('cartesian'))
    c_cart = cross / cross.norm()
    _, lat, lon = cartesian_to_spherical(c_cart.x, c_cart.y, c_cart.z)
    return SkyCoord(lon, lat, frame=c1.frame)


def cross_product_sum_skycoord2skycoord(coos):
    """
    Cross product sum of vertices, assuming vertices are in CW order.

    Use to determine minimum distance between all vertices, as a
    centroid definition.

    Parameters
    ----------
    coos : `~astropy.coordinates.SkyCoord`
        The vertices as a SkyCoord.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The minimum distance centroid cross product as a SkyCoord,
        with a spherical representation.
    """
    # verts are in CW order:
    verts_cart = coos.frame.represent_as('cartesian')
    crosssum = None

    for i in range(len(verts_cart)):
        if crosssum is None:
            crosssum = verts_cart[i - 1].cross(verts_cart[i])
        else:
            crosssum += verts_cart[i - 1].cross(verts_cart[i])

    # Normalize sum of cross products to get cartesian representation of
    # centroid === minimum distance to other points location
    c_cart = crosssum / crosssum.norm()
    _, lat, lon = cartesian_to_spherical(c_cart.x, c_cart.y, c_cart.z)

    # Ensure internal data representation format is consistent
    # with input coordinate convention:
    unit = coos[0].represent_as('spherical').lon.unit
    return SkyCoord(lon.to(unit), lat.to(unit), frame=coos.frame)


def bounding_lonlat_poles_processing(region, lons_arr, lats_arr, inner_region=None):
    """
    Check if region covers either pole & modify latitude bounds
    accordingly.

    Parameters
    ----------
    region : `~regions.SphericalSkyRegion`
        The SphericalSkyRegion region.

    lons_arr : list-like [`~astropy.coordinates.Longitude`]
        Initial estimate of longitude bounds.

    lats_arr : list-like [`~astropy.coordinates.Latitude`]
        Initial estimate of latitude bounds.

    inner_region : `~regions.SphericalSkyRegion`, optional
        The inner region, for an annulus (which modifies the pole logic).
        Default is None (for a simply-connected region).

    Returns
    -------
    lons_arr : list-like [`~astropy.coordinates.Longitude`] or None
        Corrected longitude bounds. Is None if the region goes over
        the pole (as the region contains all longitudes).

    lats_arr : list-like [`~astropy.coordinates.Latitude`]
        Corrected latitude bounds.
    """
    # Check if shape covers either pole & modify lats arr accordingly:
    poles = SkyCoord([0, 0], [-90, 90], unit=u.deg, frame=region.frame)

    # ------------------------------------
    # Simply connected region
    if inner_region is None:
        pole_contains = region.contains(poles)
        if np.any(pole_contains):
            lons_arr = None
            # S pole:
            if pole_contains[0]:
                lats_arr[0] = -90 * u.deg
            # N pole
            if pole_contains[1]:
                lats_arr[1] = 90 * u.deg

        if lons_arr is not None:
            return Longitude(lons_arr), Latitude(lats_arr)

        return None, Latitude(lats_arr)

    # ------------------------------------
    # Inner region set: annulus logic:
    pole_contains = inner_region.contains(poles)
    if np.any(pole_contains):
        lats_raw_inner = get_circle_latitude_tangent_limits(inner_region.center,
                                                            inner_region.radius)

        # S pole:
        if pole_contains[0]:
            # Change lower lats bound to be the minimum of the
            # inner boundary itself
            lats_arr[0] = lats_raw_inner[0]
        # N pole
        if pole_contains[1]:
            # Change upper lats bound to be the maximum of the
            # inner boundary itself
            lats_arr[1] = lats_raw_inner[1]

    if lons_arr is not None:
        return Longitude(lons_arr), Latitude(lats_arr)

    return None, Latitude(lats_arr)


def _get_circle_latitude_tangent_points(center, radius):
    # Get the points on an arbitrary circle that
    # intersect the tangent latitude circle limits.

    # This includes lon values, with +180deg folding
    # for over the pole cases

    crepr = center.represent_as('spherical')

    lons = [crepr.lon, crepr.lon]
    lats = [crepr.lat - radius, crepr.lat + radius]

    # Fold if over poles:
    if lats[0] < -90 * u.deg:
        lats[0] = -180 * u.deg - lats[0]
        lons[0] = lons[0] + 180 * u.deg

    if lats[1] > 90 * u.deg:
        lats[1] = 180 * u.deg - lats[1]
        lons[1] = lons[1] + 180 * u.deg

    return SkyCoord(lons, lats, frame=center.frame)


def _get_circle_longitude_tangent_points(center, radius):
    # Get the points on an arbitrary circle that
    # intersect the tangent longitude circle limits.

    # This includes lat values = 0deg
    # (all longitude circles have centers on equator)

    crepr = center.represent_as('spherical')

    lats = [0 * u.deg, 0 * u.deg]

    lats_ref = np.array(
        [(crepr.lat - radius).to(u.deg).value, (crepr.lat + radius).to(u.deg).value]
    )
    lons = []

    if np.any(np.abs(lats_ref) > 90):
        return None

    if np.any(np.abs(lats_ref) == 90):
        return SkyCoord(
            [crepr.lon - 90 * u.deg, crepr.lon + 90 * u.deg], lats, frame=center.frame
        )

    for sgn in [-1, 1]:
        # Do 1 then -1, because of lon increasing to east
        lon_gc = crepr.lon - sgn * (
            np.arccos(
                np.sin(radius.to(u.radian).value) / np.cos(crepr.lat.to(u.radian).value)
            )
            * u.radian
        ).to(u.deg)

        lons.append(lon_gc + sgn * 90 * u.deg)

    return SkyCoord(lons, lats, frame=center.frame)


def get_circle_latitude_tangent_limits(center, radius):
    """
    For an arbitrary spherical circle with center (lon0, lat0) and
    radius R0, get the tangent latitude circle limits (encoded at
    latitude values, not radii of the latitude circles).

    Use these to determine the latitude bounding limits from those tangents.
    (the points where the latitude circles intersect this circle)
    (Tangent circles have Rlat = 90 - (lat0 +- R0),
    equivalent to latitude limits lat0 +- R0)

    Ignores any issues with "over the pole" bounds --
    this only computes the values of the tangent circles.
    Other processing will handle over the pole bounds logic.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The circle center as a SkyCoord.

    radius : `~astropy.coordinates.Angle` or `~astropy.Quantity`
        The circle radius (as an |Angle| or |Quantity| with angular units).

    Returns
    -------
    latitude_limits : `~astropy.coordinates.Latitude`
        Length two |Latitude| with the latitudes of the tangent circles.
    """
    tan_lat_pts = _get_circle_latitude_tangent_points(center, radius)
    tan_lat_pts = tan_lat_pts.represent_as('spherical')

    return Latitude(tan_lat_pts.lat).to(u.deg)


def get_circle_longitude_tangent_limits(center, radius):
    """
    For an arbitrary spherical circle with center (lon0, lat0) and
    radius R0, get the longitudes of the centers of the tangent
    "longitude great circles" (as all longitude circles have
    latitude=0deg).

    Use these to determine the longitude bounds for this circle (the
    points where the longitude GCs intersect the circle).

    If | lat0 +- R0 | > 90deg: lon_bounds = None (Crosses pole, so no
    tangent longitude great circle exists -- all longitude lines [great
    circles] cross this circle.)

    If | lat0 +- R0 | = 90deg: lon_bounds = [lon0-90,lon0+90] (Touches
    either pole, so spans full 180 of longitudes centered on lon0)

    If | lat0 +- R0 | < 90deg:     lon0 +- arccos[ sin(R0) / cos(lat0) ]

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The circle center as a SkyCoord.

    radius : `~astropy.coordinates.Angle` or `~astropy.Quantity`
        The circle radius (as an |Angle| or |Quantity| with angular units).

    Returns
    -------
    longitude_limits : `~astropy.coordinates.Longitude`
        Length two |Longitude| with the longitudes of the centers
        of the tangent longitude great circles.
    """
    tan_lon_pts = _get_circle_longitude_tangent_points(center, radius)

    if tan_lon_pts is None:
        return None

    tan_lon_pts = tan_lon_pts.represent_as('spherical')
    return Longitude(tan_lon_pts.lon).to(u.deg)


def _add_tan_pts_if_in_pa_range(
    coord_list,
    tan_pts,
    gc,
    wrap_ang,
    pas_verts_wrap,
    coord=None,
):
    pas_tan_pts = gc.center.position_angle(tan_pts).to(u.deg)

    # CHECK RANGES:
    # To handle possible cases of lon values "wrapping around" across
    # the standard 360->0 wrap, wrap lon values + coord values
    # around the first entry angle
    pas_tan_pts_wrap = pas_tan_pts.wrap_at(wrap_ang)

    in_range = (pas_tan_pts_wrap >= pas_verts_wrap[0]) & (
        pas_tan_pts_wrap <= pas_verts_wrap[1]
    )

    if np.any(in_range):
        tptrepr = tan_pts.represent_as('spherical')
        coord_list = np.append(coord_list, getattr(tptrepr, coord)[in_range])

    return coord_list


def _validate_vertices_ordering(verts, gc, centroid):
    pas_verts = gc.center.position_angle(verts).to(u.deg)

    # Check ordering of vertices pas:
    pa_centroid = gc.center.position_angle(centroid).to(u.deg)
    wrap_ang = pas_verts[0]
    pas_verts_wrap = pas_verts.wrap_at(wrap_ang)
    pa_centroid_wrap = pa_centroid.wrap_at(wrap_ang)
    in_range_lon = (pa_centroid_wrap >= pas_verts_wrap[0]) & (
        pa_centroid_wrap <= pas_verts_wrap[1]
    )
    if not in_range_lon:
        # If not in range, invert the vertices PA ordering:
        pas_verts = pas_verts[::-1]
        wrap_ang = pas_verts[0]
        pas_verts_wrap = pas_verts.wrap_at(wrap_ang)

    return pas_verts_wrap, wrap_ang


def _validate_lon_bounds_ordering(lons_arr, centroid):
    # Invert longitude order if centroid is outside of range:
    wrap_ang = lons_arr[0]
    centroid_lon_wrap = centroid.represent_as('spherical').lon.wrap_at(wrap_ang)
    in_range_lon = (centroid_lon_wrap >= lons_arr[0].wrap_at(wrap_ang)) & (
        centroid_lon_wrap <= lons_arr[1].wrap_at(wrap_ang)
    )
    if not in_range_lon:
        lons_arr = lons_arr[::-1]

    return lons_arr


def get_edge_raw_lonlat_bounds_circ_edges(vertices, centroid, gcs):
    """
    Get the raw longitude / latitude bounds from the circle edges of
    spherical sky region.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices as a SkyCoord.

    centroid : `~astropy.coordinates.SkyCoord`
        The polygon centroid as a SkyCoord.

    gcs : `~regions.CircleSphericalSkyRegion`
        The circle boundaries as CircleSphericalSkyRegion instances.

    Returns
    -------
    longitude_limits, latitude_limits: `~astropy.coordinates.Longitude`
        Length two |Longitude| and |Latitude| with the computed
        longitude/latitude bounds from the polygon edges.
    """
    # Consider lon/lat of vertices: may produce min/max bounds:
    vrepr = vertices.represent_as('spherical')

    # lons_list = vrepr.lon
    # lats_list = vrepr.lat

    # Special handling:
    # Exclude vertices from longitude bounds if any is on a pole
    lons_list = []
    lats_list = []
    for v in vrepr:
        if np.abs(v.lat.to(u.deg).deg) < 90:
            lons_list.append(v.lon)
            lats_list.append(v.lat)
    lons_list = Longitude(lons_list, unit=u.radian)
    lats_list = Latitude(lats_list, unit=u.radian)

    # Need to also check for "out bulging" from edges,
    # as far as latitude/lon bounds:
    # eg, 2 vertices at ~60deg: the gc arc goes ~closer to the pole;
    # a circle centered close to the pole but not extending over it:
    # WIDE lon bounds

    # Requires checking if the "bounding" points are *ON* the polygons.

    for i, gc in enumerate(gcs):
        # PAs from gc center to vertices:
        verts = concatenate([vertices[i - 1], vertices[i]])

        pas_verts_wrap, wrap_ang = _validate_vertices_ordering(verts, gc, centroid)

        # --------------------------------------------------------
        # Latitude tangent points from bound circle as len 2 SkyCoord:
        # Only add to the list if the tangent point is located along this edge

        tan_lat_pts = _get_circle_latitude_tangent_points(gc.center, gc.radius)

        lats_list = _add_tan_pts_if_in_pa_range(
            lats_list, tan_lat_pts, gc, wrap_ang, pas_verts_wrap, coord='lat'
        )

        # --------------------------------------------------------
        # Longitude tangent points from bound circle as len 2 SkyCoord:
        # Only add to the list if the tangent point is located along this edge

        tan_lon_pts = _get_circle_longitude_tangent_points(gc.center, gc.radius)
        if tan_lon_pts is not None:
            lons_list = _add_tan_pts_if_in_pa_range(
                lons_list, tan_lon_pts, gc, wrap_ang, pas_verts_wrap, coord='lon'
            )

    lons_arr = [lons_list.min(), lons_list.max()]
    lats_arr = [lats_list.min(), lats_list.max()]

    # --------------------------------------------------------
    # Invert longitude order if centroid is outside of range:
    lons_arr = _validate_lon_bounds_ordering(lons_arr, centroid)

    return Longitude(lons_arr).to(u.deg), Latitude(lats_arr).to(u.deg)


def _discretize_edge_boundary(vertices, circ, centroid, n_points):
    # Discretize an edge boundary defined by a circle, geodetic or not:
    # either great circle arc, or a span of a non-great circle
    # (e.g., constant lat edges of RangeSphericalSkyRegion)

    # For every edge boundary: determine range of PAs spanned by lines
    # connecting circle center to the two vertices bounding that edge:
    pas_verts = circ.center.position_angle(vertices).to(u.deg)

    pas_verts_wrap, wrap_ang = _validate_vertices_ordering(vertices, circ, centroid)

    # Need to wrap angles to calculate span: wrap at lower value:
    pas_verts_wrap = pas_verts.wrap_at(wrap_ang)
    pa_span = pas_verts_wrap[1] - pas_verts_wrap[0]

    # Sample angle range over pa_span with Npoints, and offset
    # to start at pas_verts[0]
    theta = np.linspace(0, 1, num=n_points, endpoint=False) * pa_span + pas_verts[0]

    # Calculate directional offsets to get boundary discretization,
    # with vertices as SkyCoords
    bound_verts = circ.center.directional_offset_by(theta, circ.radius)

    return bound_verts


def discretize_all_edge_boundaries(vertices, circs, centroid, n_points):
    """
    Discretize all edge boundaries for spherical sky regions.

    Parameters
    ----------
    vertices : `~astropy.coordinates.SkyCoord`
        The vertices of the spherical sky region.

    circs : `~regions.CircleSphericalSkyRegion`
        The circle boundaries as CircleSphericalSkyRegion instances.

    centroid : `~astropy.coordinates.SkyCoord`
        The region centroid as a SkyCoord.

    n_points : int
        The number of points for discretization along each edge.

    Returns
    -------
    all_edge_bound_verts : `~astropy.coordinates.SkyCoord`
        The discretized boundary edge vertices.
    """
    # Iterate over full set of vertices & boundary circles for
    # a region (polygon or range)

    all_edge_bound_verts = None
    for i, circ in enumerate(circs):
        # PAs from gc center to vertices:
        verts = concatenate([vertices[i - 1], vertices[i]])

        bound_verts = _discretize_edge_boundary(verts, circ, centroid, n_points)

        if all_edge_bound_verts is None:
            all_edge_bound_verts = bound_verts
        else:
            all_edge_bound_verts = concatenate(
                [all_edge_bound_verts.copy(), bound_verts]
            )
            # For some reason concatenate is adding distances,
            # so use remove those by running from
            # UnitSpherical->SphericalRepresentation...
            all_edge_bound_verts = SkyCoord(
                SkyCoord(
                    all_edge_bound_verts,
                    representation_type=UnitSphericalRepresentation,
                ),
                representation_type=SphericalRepresentation,
            )

    return all_edge_bound_verts
