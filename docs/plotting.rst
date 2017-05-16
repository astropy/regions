.. _gs-mpl:

Plotting regions with Matplotlib
================================

Some `~regions.PixelRegion` objects have an ``as_patch()`` method that returns an
equivalent `matplotlib.patches` object. For example :meth:`regions.CirclePixelRegion.as_patch`
returns a `matplotlib.patches.Circle` object.

To draw a matplotlib patch object, add it to an `matplotlib.axes.Axes` object.

.. plot::
   :include-source:

    from regions import PixCoord, CirclePixelRegion
    import matplotlib.pyplot as plt

    region = CirclePixelRegion(PixCoord(x=0.3, y=0.42), radius=0.5)
    patch = region.as_patch()

    axes = plt.gca()
    axes.set_aspect('equal')
    axes.add_patch(patch)

    plt.show()


The :meth:`~regions.PixelRegion.plot` convenience method just does these two
steps at once (creating a matplotlib patch and adding it to an axis),
and does call ``plt.gca()`` if no axis is passed in.

Note that not all pixel regions have ``as_patch()`` methods, e.g.
the `~regions.PointPixelRegion` or compound regions don't because there's
no equivalent matplotlib object.

Here's a full example how to plot a `~regions.CirclePixelRegion` on an image.

.. plot:: plot_example_pix.py
   :include-source:

The `~regions.RectanglePixelRegion` and `~regions.EllipsePixelRegion` docstrings also
contain plot examples.

`~regions.SkyRegion` objects currently don't have an ``as_patch()`` or ``plot()``
method. To plot them, convert them to a pixel region first:

.. code-block:: python

    sky_region = <...>
    pixel_region = sky_region.to_pixel(wcs, mode, tolerance)
    pixel_region.plot(**kwargs)  # plot options passed to matplotlib

We do plan to add extensive documentation on sky region plotting, or to
add methods on sky region to do it directly in the future
(see https://github.com/astropy/regions/issues/76 ),
after the polygon region classes are developed.

An example of how to plot sky regions on a sky image is shown above.
