Plotting Regions with Matplotlib
================================

Some `~regions.PixelRegion` objects have an ``as_artist()``
method that returns an equivalent `matplotlib.patches` object.
For example :meth:`regions.CirclePixelRegion.as_artist` returns a
`matplotlib.patches.Circle` object.

To draw a matplotlib patch object, add it to an `matplotlib.axes.Axes`
object.

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

The :meth:`regions.PixelRegion.plot` method is a convenience method that
combines these two steps (creating a matplotlib patch artist and adding
it to an axis). If no axis is passed then it calls ``plt.gca()``.

You can shift the origin of the region while plotting by supplying the
``origin`` pixel coordinates to either :meth:`~regions.PixelRegion.plot`
or :meth:`~regions.PixelRegion.as_artist`. The
:meth:`~regions.PixelRegion.plot` method also takes any keyword argument
that the `~matplotlib.patches.Patch` object accepts for those regions
represented as patches, or arguments `~matplotlib.lines.Line2D` accepts
for point regions. For example:

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from regions import PixCoord, CirclePixelRegion

    fig, ax = plt.subplots()
    region = CirclePixelRegion(center=PixCoord(x=7, y=5), radius=3)

    data = np.arange(10 * 15).reshape((10, 15))
    ax.imshow(data, cmap='gray', interpolation='nearest', origin='lower')
    region.plot(ax=ax, color='red', lw=2.0)

The documentation for `~regions.RectanglePixelRegion` and
`~regions.EllipsePixelRegion` also shows plotting examples.


Plotting Sky Regions
--------------------

Note that `~regions.SkyRegion` objects do not have an ``as_artist()`` or
``plot()`` method. To plot a `~regions.SkyRegion` object, you will need
to convert it to a pixel region (using a WCS object):

.. doctest-skip::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion
    >>> sky_center = SkyCoord(42, 43, unit='deg')
    >>> sky_radius = Angle(25, 'deg')
    >>> sky_region = CircleSkyRegion(sky_center, sky_radius)
    >>> pixel_region = sky_region.to_pixel(wcs)
    >>> pixel_region.plot()
