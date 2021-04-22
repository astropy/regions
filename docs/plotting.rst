.. _gs-mpl:

Plotting regions with Matplotlib
================================

Some `~regions.PixelRegion` objects have an ``as_artist()`` method that returns an
equivalent `matplotlib.patches` object. For example :meth:`regions.CirclePixelRegion.as_artist`
returns a `matplotlib.patches.Circle` object.

To draw a matplotlib patch object, add it to an `matplotlib.axes.Axes` object.

.. plot::
   :include-source:

    from regions import PixCoord, CirclePixelRegion
    import matplotlib.pyplot as plt

    region = CirclePixelRegion(PixCoord(x=0.3, y=0.42), radius=0.5)
    artist = region.as_artist()

    axes = plt.gca()
    axes.set_aspect('equal')
    axes.add_artist(artist)
    axes.set_xlim([-0.5, 1])
    axes.set_ylim([-0.5, 1])


The :meth:`~regions.PixelRegion.plot`, a convenience method just does these two
steps at once (creating a matplotlib patch artist and adding it to an axis),
and calls ``plt.gca()`` if no axis is passed in.

You can shift the origin of the region very conveniently while plotting by simply
supplying the ``origin`` pixel coordinates to :meth:`~regions.PixelRegion.plot`
and :meth:`~regions.PixelRegion.as_artist`. The ``**kwargs`` argument takes any
keyword argument that the `~matplotlib.patches.Patch` object accepts for those
regions represented as patches, or arguments `~matplotlib.lines.Line2D` accepts
for point regions. For example:

.. plot::
   :include-source:

    from regions import PixCoord, BoundingBox
    import matplotlib.pyplot as plt

    bbox = BoundingBox(ixmin=-1, ixmax=1, iymin=-2, iymax=2)
    # shifting the origin to (1, 1) pixel position
    ax = bbox.plot(origin=(1, 1), edgecolor='yellow', facecolor='red', fill=True)
    ax.set_xlim([-4, 2])
    ax.set_ylim([-4, 2])


Here's a full example how to plot a `~regions.CirclePixelRegion` on an image.

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from regions import PixCoord, CirclePixelRegion

    fig, ax = plt.subplots()
    region = CirclePixelRegion(center=PixCoord(x=3, y=5), radius=3)

    data = np.arange(10 * 15).reshape((10, 15))
    ax.imshow(data, cmap='gray', interpolation='nearest', origin='lower')
    region.plot(ax=ax, color='red')


The `~regions.RectanglePixelRegion` and `~regions.EllipsePixelRegion` docstrings also
contain plot examples.

`~regions.SkyRegion` objects currently don't have an ``as_artist()`` or ``plot()``
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
