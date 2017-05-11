.. _gs-shapely:

Exporting regions to shapely objects
====================================

The `Shapely <http://toblerity.org/shapely/manual.html>`__ Python package is a
generic package for the manipulation and analysis of geometric objects in the
Cartesian plane. Concerning regions in the cartesian plane, it is more
feature-complete, powerful and optimized than this ``regions`` package. It
doesn't do everything astronomers need though, e.g. no sky regions, no use of
Astropy classes like `~astropy.coordinates.SkyCoord` or
`~astropy.coordinates.Angle`, and no region serialisation with the formats
astronomers use (such as e.g. ds9 regions).

`~regions.PixelRegion` classes provide a method :meth:`~regions.PixelRegion.to_shapely` that allows creation
of Shapely shape objects. At the moment there is no ``from_shapely`` method to convert Shapely objects
back to ``regions`` objects. The future of the use of Shapely in ``regions`` is currently unclear, some options are:

1. Add ``from_shapely`` and use it to implement e.g. `~regions.PolygonPixelRegion` operations
   can "discretization" of other shapes to polygons.
   This would make Shapely a required dependency to work with polygons.
2. Keep ``to_shapely`` for the (small?) fraction of users that want to do this,
   but don't expand or use it inside the ``regions`` package  to avoid the extra heavy dependency.
3. Remove the use of Shapely completely from the API unless good use cases demonstrating a need come up.

Here's an example how to create a Shapely object and do something that's not implemented in ``regions``,
namely to buffer a rectangle, resulting in a polygon.

.. code-block:: python

    # TODO: RectanglePixelRegion isn't implemented yet, this doesn't work yet.
    from regions import RectanglePixelRegion
    region = RectanglePixelRegion(center=(3, 2), width=2, height=1)
    shape = region.to_shapely()
    # `shape` is a `shapely.geometry.polygon.Polygon` object
    shape2 = shape.buffer(distance=3)
    # `shape2` is a polygon that's buffered by 3 pixels compared to `shape`
