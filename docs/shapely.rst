Converting Regions to Shapely Objects
=====================================

The `Shapely <https://shapely.readthedocs.io/en/latest/>`__ Python
package is a generic package for the manipulation and analysis of
geometric objects in the Cartesian plane. For regions in the cartesian
plane, it is more feature-complete, powerful, and optimized than this
package.

The use of Shapely or other Python regions packages that come from the
geospatial domain in Astronomy is rare. However, if you have a complex
pixel region analysis task, you can consider using Shapely. Either use
it directly, by defining Shapely regions via Python code or one of the
serialization formats they support, or by writing some Python code to
convert `~regions.Region` objects to Shapely objects.

Here we give an example of converting a circle to a Shapely object
and then polygonize it. One nice feature of Shapely is that it can
polygonize all shapes and perform fast polygon-based computations like
intersection and union.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from regions import PixCoord, CirclePixelRegion
    from shapely.geometry import Point

    # Make an example region
    region = CirclePixelRegion(center=PixCoord(3, 2), radius=2)

    # Convert to Shapely
    point = Point(region.center.x, region.center.y)
    circle = point.buffer(region.radius)

    # Actually, this is a polygon approximation of the circle
    print(circle)

    # Plot the result
    x, y = circle.exterior.xy
    ax = plt.subplot(1, 1, 1)
    ax.set_aspect('equal')
    ax.plot(x, y, 'g-')
