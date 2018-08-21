.. _gs-shapely:

Converting regions to shapely objects
=====================================

The `Shapely <http://toblerity.org/shapely/manual.html>`__ Python package is a
generic package for the manipulation and analysis of geometric objects in the
Cartesian plane. Concerning regions in the cartesian plane, it is more
feature-complete, powerful and optimized than this ``regions`` package.

The use of Shapely or other Python regions packages that come from the geospatial domain
in Astronomy is rare. However, if you have a complex pixel region analysis task,
you can consider using Shapely. Either use it directly, by defining Shapely regions
via Python code or one of the serialisation formats they support, or by writing
some Python code to convert ``astropy-regions`` objects to Shapely objects.

Here we give one example how to do this: convert a circle to a Shapely object
and polygonise it. That's one nice feature of Shapely, it can polygonise all shapes
and do fast polygon-based computations like intersection and union. If you need to
do this, that's a good reason to use Shapely.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from regions import PixCoord, CirclePixelRegion

    # Make an example region
    region = CirclePixelRegion(center=PixCoord(3, 2), radius=2)

    # Convert to Shapely
    from shapely.geometry import Point
    point = Point(region.center.x, region.center.y)
    circle = point.buffer(region.radius)

    # Actually, this is a polygon approximation of the circle!
    print(circle)

    # Plot the result
    x, y = circle.exterior.xy
    ax = plt.subplot(1, 1, 1)
    ax.plot(x, y, 'g-')
    plt.show()
