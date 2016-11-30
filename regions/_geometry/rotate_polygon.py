# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy import units as u
from astropy.coordinates.representation import UnitSphericalRepresentation

try:
    from astropy.coordinates.matrix_utilities import rotation_matrix as rotation_array
    def rotation_matrix(*args, **kwargs):
        return np.matrix(rotation_array(*args, **kwargs))
except ImportError:  # Astropy < 1.3
    from astropy.coordinates.angles import rotation_matrix

__all__ = ['rotate_polygon']


def rotate_polygon(lon, lat, lon0, lat0):
    """
    Given a polygon with vertices defined by (lon, lat), rotate the polygon
    such that the North pole of the spherical coordinates is now at (lon0,
    lat0). Therefore, to end up with a polygon centered on (lon0, lat0), the
    polygon should initially be drawn around the North pole.
    """

    # Create a representation object
    polygon = UnitSphericalRepresentation(lon=lon, lat=lat)

    # Determine rotation matrix to make it so that the circle is centered
    # on the correct longitude/latitude.
    m1 = rotation_matrix(-(0.5 * np.pi * u.radian - lat0), axis='y')
    m2 = rotation_matrix(-lon0, axis='z')
    transform_matrix = m2 * m1

    # Apply 3D rotation
    polygon = polygon.to_cartesian()
    polygon = polygon.transform(transform_matrix)
    polygon = UnitSphericalRepresentation.from_cartesian(polygon)

    return polygon.lon, polygon.lat
