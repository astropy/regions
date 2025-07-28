# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides spherical sky region calculation helper tools.
"""

import astropy.units as u
import numpy as np
from astropy.coordinates import (SkyCoord, SphericalRepresentation,
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

    return SkyCoord(lon, lat, frame=coos.frame)


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
