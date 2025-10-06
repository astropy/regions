
.. _spherical_bounding_regions:

Bounding Regions for On-Sky Spatial Searches
============================================

The `~regions.SphericalSkyRegion` classes also provide
`~regions.SphericalSkyRegion.bounding_circle` and
`~regions.SphericalSkyRegion.bounding_lonlat` properties,
yielding the spherical circle and longitude and latitude ranges
enclosing the region, respectively.

This bounding information can be helpful when querying databases
to select objects with tools such as
`Astroquery <https://astroquery.readthedocs.io/en/latest/>`_
or `PyVO's TAP service <https://pyvo.readthedocs.io/en/latest/dal/index.html#pyvo-tap>`_.

For many regions of interest, the spherical region properties can
directly provide the query search boundary
(e.g., the center and radius for a circle for a cone search,
or the vertices for a polygon region search).
However, this depends on what query shapes are supported by
the database.

When a particular region shape is not supported by a
database, these bounding properties can be used to define a
"padded" search, and the results can then be downselected to only
those within the region using the `~regions.SphericalSkyRegion.contains()`.


(In particular, the bounding properties may be helpful when a
coordinate frame transformation is necessary to match the database
coordinate system; see :ref:`spherical_frame_transform`.)


As an example, let's create two spherical regions, a circle and range,
and determine their bounding longitude/latitude spans and circle:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import CircleSphericalSkyRegion, RangeSphericalSkyRegion
    >>> sph_circ = CircleSphericalSkyRegion(SkyCoord(100,-40,unit=u.deg,
    ...                                              frame="galactic"),
    ...                                     30*u.deg)
    >>> print(sph_circ.bounding_lonlat)
    (<Longitude [ 59.25424338, 140.74575662] deg>, <Latitude [-70., -10.] deg>)
    >>> sph_range = RangeSphericalSkyRegion(frame="galactic",
    ...                                     longitude_range=[-45,45]*u.deg,
    ...                                     latitude_range=[0,45]*u.deg)
    >>> print(sph_range.bounding_circle)
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (0., 16.32494994)>
    radius: 47.2657903991002 deg

(The bounding circle of a circle is always itself, and in the original
coordinate frame, the bounding longitude/latitude of a
spherical range is equivalent to the region.)

Let's now transform these regions into the ICRS frame,
and obtain the new bounds, including the new frame longitude/latitude
bounds on the transformed Range region:


.. code-block:: python

    >>> sph_circ_transf = sph_circ.transform_to("icrs")
    >>> print(sph_circ_transf.bounding_lonlat)
    (<Longitude [322.34509896,  26.43985125] deg>, <Latitude [-10.44033547,  49.55966453] deg>)
    >>> sph_range_transf = sph_range.transform_to("icrs")
    >>> print(sph_range_transf.bounding_lonlat)
    (<Longitude [201.75437889, 288.42757587] deg>, <Latitude [-60.49575721,  27.00079109] deg>)
    >>> print(sph_range_transf.bounding_circle)
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (251.64674219, -19.64298539)>
    radius: 47.26579039910021 deg


These bounding circles and longitude/latitude are visualized below for both
the original and transformed coordinate frames.


.. plot::
   :include-source: false

    from astropy.coordinates import SkyCoord, Angle
    from astropy.visualization.wcsaxes.frame import EllipticalFrame
    from astropy import units as u

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from regions import (CircleSphericalSkyRegion,
                         RangeSphericalSkyRegion,
                         make_example_dataset)


    dataset = make_example_dataset(data='simulated')
    wcs = dataset.wcs

    dataset_icrs = make_example_dataset(data='simulated', config={'ctype':
                                                             ('RA---AIT',
                                                              'DEC--AIT')})
    wcs_icrs = dataset_icrs.wcs

    sph_circ = CircleSphericalSkyRegion(SkyCoord(100,-40,
                                                 unit=u.deg,
                                                 frame="galactic"),
                                        30*u.deg)
    sph_range = RangeSphericalSkyRegion(frame="galactic",
                                        longitude_range=[-45,45]*u.deg,
                                        latitude_range=[0,45]*u.deg)
    sph_circ_transf = sph_circ.transform_to("icrs")
    sph_range_transf = sph_range.transform_to("icrs")


    fig = plt.figure()
    fig.set_size_inches(7,7)

    axes = []
    axes.append(fig.add_axes([0.15, 0.575, 0.8, 0.425],
                             projection=wcs,
                             frame_class=EllipticalFrame))
    axes.append(fig.add_axes([0.15, 0.05, 0.8, 0.425],
                             projection=wcs_icrs,
                             frame_class=EllipticalFrame))

    ax = axes[0]
    ax.coords.grid(color='gray')
    ax.set_xlabel(r"Galactic $\ell$", labelpad=10)
    ax.set_ylabel(r"Galactic $b$")
    ax.set_title("Galactic coordinates", pad=5)

    patch = sph_circ.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, color='tab:blue', lw=3)

    sph_range.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":250}
    ).plot(ax=ax, color='tab:red', lw=3)

    bound_color = 'tab:blue'
    bound_lw = 0.75
    bound_zord = 2
    bound_lon, bound_lat = sph_circ.bounding_lonlat
    for i in range(2):
        # Need to manually "super sample" to get correct distortion.
        # Unfortunately axv/hline works for just "aitoff" projection,
        # not doing a WCS it seems.
        npt = 250
        xarr = np.repeat(bound_lon[i].degree, npt)
        yarr = np.linspace(-90,90,num=npt,endpoint=True)
        l1 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('galactic'))
        xarr = np.linspace(-180,180,num=npt,endpoint=True)
        yarr = np.repeat(bound_lat[i].degree, npt)
        l2 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('galactic'))
        ax.add_artist(l1)
        ax.add_artist(l2)


    bound_color = 'tab:red'
    sph_range.bounding_circle.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, color='tab:red', lw=bound_lw, ls='--', zorder=bound_zord)

    patch.set_clip_path(ax.coords.frame.patch)

    ax.set_xlim(20,340)
    ax.set_ylim(10,170)

    ax = axes[1]
    ax.coords.grid(color='gray')
    ax.set_xlabel("RA", labelpad=10)
    ax.set_ylabel("Dec")
    ax.set_title("ICRS coordinates", pad=5)

    patch = sph_circ_transf.to_pixel(
       wcs=wcs_icrs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, color='tab:blue', lw=3)

    sph_range_transf.to_pixel(
       wcs=wcs_icrs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":250}
    ).plot(ax=ax, color='tab:red', lw=3)


    bound_color = 'tab:blue'
    bound_lw = 0.75
    bound_zord = 2
    bound_lon, bound_lat = sph_circ_transf.bounding_lonlat
    for i in range(2):
        # Need to manually "super sample" to get correct distortion.
        # Unfortunately axv/hline works for just "aitoff" projection,
        # not doing a WCS it seems.
        npt = 250
        xarr = np.repeat(bound_lon[i].degree, npt)
        yarr = np.linspace(-90,90,num=npt,endpoint=True)
        l1 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('icrs'))
        xarr = np.linspace(-180,180,num=npt,endpoint=True)
        yarr = np.repeat(bound_lat[i].degree, npt)
        l2 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('icrs'))
        ax.add_artist(l1)
        ax.add_artist(l2)


    bound_color = 'tab:red'
    bound_lon, bound_lat = sph_range_transf.bounding_lonlat
    for i in range(2):
        # Need to manually "super sample" to get correct distortion.
        # Unfortunately axv/hline works for just "aitoff" projection,
        # not doing a WCS it seems.
        npt = 250
        xarr = np.repeat(bound_lon[i].degree, npt)
        yarr = np.linspace(-90,90,num=npt,endpoint=True)
        l1 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('icrs'))
        xarr = np.linspace(-180,180,num=npt,endpoint=True)
        yarr = np.repeat(bound_lat[i].degree, npt)
        l2 = Line2D(xarr, yarr, ls='--', color=bound_color,
                    lw=bound_lw, zorder=bound_zord,
                    transform=ax.get_transform('icrs'))
        ax.add_artist(l1)
        ax.add_artist(l2)
    sph_range_transf.bounding_circle.to_pixel(
       wcs=wcs_icrs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, color='tab:red', lw=bound_lw, ls='--', zorder=bound_zord)


    patch.set_clip_path(ax.coords.frame.patch)

    ax.set_xlim(20,340)
    ax.set_ylim(10,170)
    ax.coords[0].set_format_unit(u.deg)
