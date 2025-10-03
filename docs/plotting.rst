Plotting Regions with Matplotlib
================================

Plotting Pixel Regions
----------------------

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

The default keyword arguments for the matplotlib artist depend on the
value of the ``default_style`` keyword in the `~regions.RegionVisual`
dictionary. This keyword is currently set (to a value of 'ds9') only
when reading from DS9 region files. If this keyword is not set or set
to 'mpl' or `None`, then the matplotlib defaults will be used, with the
exception that fill is turned off for `~matplotlib.patches.Patch` and
`~matplotlib.lines.Line2D` artists.

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


Plotting Spherical Sky Regions
------------------------------

Similarlly, `~regions.SphericalSkyRegion` objects do not have an
``as_artist()`` or ``plot()`` method. To plot a `~regions.SphericalSkyRegion`
object, you will need to convert it to a pixel region (using a WCS object).
Boundary distortions effects can also be included in this conversion (by
setting the ``include_boundary_distortions`` keyword), to capture
the effects of projecting from spherical to a planar geometry.
(See the second example in :ref:`index_examples`.)

It is also possible to use the coordinates of a discretized
`~regions.SphericalSkyRegion` to show the region's boundary in a figure.

.. plot::
   :include-source:

    from astropy.coordinates import Angle, SkyCoord
    import matplotlib.pyplot as plt

    from regions import CircleSphericalSkyRegion, make_example_dataset

    dataset = make_example_dataset(data='simulated')
    wcs = dataset.wcs

    sph_sky_center = SkyCoord(42, 30, unit='deg', frame='galactic')
    sph_sky_radius = Angle(12, 'deg')
    sph_sky_region = CircleSphericalSkyRegion(sph_sky_center, sph_sky_radius)

    fig, ax = plt.subplots(figsize=(8,4),
                           subplot_kw=dict(projection=wcs))
    ax.grid(True)
    ax.set_xlabel(r"Galactic $\ell$")
    ax.set_ylabel(r"Galactic $b$")

    sph_sky_region.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, color='tab:red', lw=3)

    sph_sky_center2 = SkyCoord(42, 43, unit='deg', frame='galactic')
    sph_sky_radius2 = Angle(25, 'deg')
    sph_sky_region2 = CircleSphericalSkyRegion(sph_sky_center2, sph_sky_radius2)
    poly_sph_sky2 = sph_sky_region2.discretize_boundary(n_points=1000)
    ax.plot(
        poly_sph_sky2.vertices.l,
        poly_sph_sky2.vertices.b,
        lw=2, color="tab:blue",
        transform=ax.get_transform('galactic'),
    )
