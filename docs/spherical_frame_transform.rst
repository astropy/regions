.. include:: references.txt

.. _spherical_frame_transform:

Coordinate Frame Transformation of Spherical Regions
====================================================

The `~astropy.coordinates.SkyCoord` class provides functionality for
representing celestial coordinates and for transforming between
different coordinate frames (see the :ref:`astropy-coordinates` documentation
for full details).

Since the `~regions.SphericalSkyRegion` class represents geometric shapes on the
celestial sphere, these idealized shapes can also be represented in, and
transformed between, any spherical coordinate reference frame.
Transforming spherical regions can be useful in cases where a region of
interest is defined in one coordinate frame, but a query
searching for targets in this region needs to be specified
in second coordinate frame.

As an example, let's start by defining two sky regions in the Galactic
coordinate frame:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import CircleSphericalSkyRegion, RangeSphericalSkyRegion
    >>> sph_circ = CircleSphericalSkyRegion(SkyCoord(100,-30,unit=u.deg,
    ...                                              frame="galactic"),
    ...                                     30*u.deg)
    >>> sph_range = RangeSphericalSkyRegion(frame="galactic",
    ...                                     longitude_range=[315,45]*u.deg,
    ...                                     latitude_range=[0,45]*u.deg)
    >>> print(sph_circ)  # doctest: +FLOAT_CMP
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (100., -30.)>
    radius: 30.0 deg
    >>> print(sph_range)  # doctest: +FLOAT_CMP
    Region: RangeSphericalSkyRegion
    frame: galactic
    longitude_range: [315.  45.] deg
    latitude_range: [ 0. 45.] deg



To convert these regions to the ICRS coordinate frame, the
`~regions.SphericalSkyRegion.transform_to()` method is used.

.. code-block:: python

    >>> sph_circ_transf = sph_circ.transform_to("icrs")
    >>> sph_range_transf = sph_range.transform_to("icrs")
    >>> print(sph_circ_transf)  # doctest: +FLOAT_CMP
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (350.21026136, 28.80607705)>
    radius: 30.0 deg
    >>> print(sph_range_transf)  # doctest: +FLOAT_CMP
    Region: RangeSphericalSkyRegion
    frame: icrs
    longitude_bounds: <LuneSphericalSkyRegion(center_gc1=<SkyCoord (ICRS): (ra, dec) in deg
        (108.42757587, -10.72370325)>, center_gc2=<SkyCoord (ICRS): (ra, dec) in deg
        (37.98010947, 60.49575721)>)>
    latitude_bounds: <CircleAnnulusSphericalSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in deg
        (192.85947789, 27.12825241)>, inner_radius=45.0 deg, outer_radius=90.0 deg)>


Note that the Range region boundaries cannot be simply described
by longitude/latitude boundaries in a transformed frame, so
the underlying circular annulus and lune boundaries (capturing
the original frame latitude and longitude bounds, respectively)
are used to describe this Range regions after a transformation.

These original and transformed regions are shown below on
full-sky projections for the Galactic and ICRS coordinate
reference frames (respectively).

.. plot::
   :include-source: false

    from astropy.coordinates import SkyCoord, Angle
    from astropy.visualization.wcsaxes.frame import EllipticalFrame
    from astropy import units as u

    import matplotlib.pyplot as plt

    from regions import (CircleSphericalSkyRegion,
                         RangeSphericalSkyRegion,
                         make_example_dataset)


    dataset = make_example_dataset(data='simulated')
    wcs = dataset.wcs

    dataset_icrs = make_example_dataset(data='simulated', config={'ctype':
                                                             ('RA---AIT',
                                                              'DEC--AIT')})
    wcs_icrs = dataset_icrs.wcs

    sph_circ = CircleSphericalSkyRegion(SkyCoord(100,-30,
                                                 unit=u.deg,
                                                 frame="galactic"),
                                        30*u.deg)
    sph_range = RangeSphericalSkyRegion(frame="galactic",
                                        longitude_range=[315,45]*u.deg,
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
    ax.coords.grid(color='black')
    ax.set_xlabel(r"Galactic $\ell$", labelpad=10)
    ax.set_ylabel(r"Galactic $b$")
    ax.set_title("Galactic coordinates", pad=5)

    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='gray', ls='dotted')

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

    patch.set_clip_path(ax.coords.frame.patch)

    ax.set_xlim(20,340)
    ax.set_ylim(10,170)

    ax = axes[1]
    ax.coords.grid(color='gray', ls='dotted')
    ax.set_xlabel("RA", labelpad=10)
    ax.set_ylabel("Dec")
    ax.set_title("ICRS coordinates", pad=5)

    overlay = ax.get_coords_overlay('galactic')
    overlay.grid(color='black', ls='solid')

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

    patch.set_clip_path(ax.coords.frame.patch)

    ax.set_xlim(20,340)
    ax.set_ylim(10,170)
    ax.coords[0].set_format_unit(u.deg)
