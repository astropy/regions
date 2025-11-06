Combining Regions
=================

There are a few ways to combine any two `~regions.Region` objects
into a compound region, i.e., a `~regions.CompoundPixelRegion`,
`~regions.CompoundSkyRegion`, or `~regions.CompoundSphericalSkyRegion` object.

Let's start by defining two sky regions::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion

    >>> circle1 = CircleSkyRegion(
    ...    center=SkyCoord(1, 2, unit='deg', frame='galactic'),
    ...    radius=Angle('5 deg'))
    >>> circle2 = CircleSkyRegion(
    ...    center=SkyCoord(-4, 3, unit='deg', frame='galactic'),
    ...    radius=Angle('3 deg'))


Intersection (AND)
------------------

To create an intersection compound region, use either the ``&`` operator
or the :meth:`~regions.Region.intersection` method::

    >>> comp_region = circle1 & circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function and_>

    >>> comp_region = circle1.intersection(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function and_>


Union (OR)
----------

To create a union compound region, use either the ``|`` operator or the
:meth:`~regions.Region.union` method::

    >>> comp_region = circle1 | circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function or_>

    >>> comp_region = circle1.union(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function or_>


Symmetric Difference (XOR)
--------------------------

To create a symmetric difference compound region, use either the ``^``
operator or the :meth:`~regions.Region.symmetric_difference` method::

    >>> comp_region = circle1 ^ circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function xor>

    >>> comp_region = circle1.symmetric_difference(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function xor>


Example Illustrating Compound Regions
-------------------------------------

The following examples demonstrate how to combine planar sky regions
and spherical sky regions, with the same circle centers and radii.

.. plot::
    :include-source: false

    # planar

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.coordinates import Angle, SkyCoord

    from regions import CircleSkyRegion, make_example_dataset

    # load example dataset to get skymap
    config = dict(crval=(0, 0),
                crpix=(180, 90),
                cdelt=(-1, 1),
                shape=(180, 360))

    dataset = make_example_dataset(data='simulated', config=config)
    wcs = dataset.wcs

    # remove sources
    dataset.image.data = np.zeros_like(dataset.image.data)

    #----------------------------------------
    # define 2 sky circles
    circle1 = CircleSkyRegion(
        center=SkyCoord(20, 0, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))

    circle2 = CircleSkyRegion(
        center=SkyCoord(50, 45, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))

    # define skycoords
    lon = np.arange(-180, 181, 10)
    lat = np.arange(-90, 91, 10)
    coords = np.array(np.meshgrid(lon, lat)).T.reshape(-1, 2)
    skycoords = SkyCoord(coords, unit='deg', frame='galactic')

    #----------------------------------------
    # get events in AND and XOR
    compound_and = circle1 & circle2
    compound_xor = circle1 ^ circle2

    mask_and = compound_and.contains(skycoords, wcs)
    skycoords_and = skycoords[mask_and]
    mask_xor = compound_xor.contains(skycoords, wcs)
    skycoords_xor = skycoords[mask_xor]

    # plot
    fig = plt.figure()
    fig.set_size_inches(7,3.5)
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs, aspect='equal')

    ax.scatter(skycoords.l.value, skycoords.b.value, label='all',
            transform=ax.get_transform('galactic'))
    ax.scatter(skycoords_xor.l.value, skycoords_xor.b.value, color='orange',
            label='xor', transform=ax.get_transform('galactic'))
    ax.scatter(skycoords_and.l.value, skycoords_and.b.value, color='magenta',
            label='and', transform=ax.get_transform('galactic'))

    circle1.to_pixel(wcs=wcs).plot(ax=ax, edgecolor='green', facecolor='none',
                                alpha=0.8, lw=3)
    circle2.to_pixel(wcs=wcs).plot(ax=ax, edgecolor='red', facecolor='none',
                                alpha=0.8, lw=3)

    ax.legend(loc='lower right')

    ax.set_xlim(-0.5, dataset.config['shape'][1] - 0.5)
    ax.set_ylim(-0.5, dataset.config['shape'][0] - 0.5)
    ax.set_title("Planar SkyRegions")


.. plot::
    :include-source: false

    # spherical

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.coordinates import Angle, SkyCoord

    from regions import CircleSphericalSkyRegion, make_example_dataset

    # load example dataset to get skymap
    config = dict(crval=(0, 0),
                crpix=(180, 90),
                cdelt=(-1, 1),
                shape=(180, 360))

    dataset = make_example_dataset(data='simulated', config=config)
    wcs = dataset.wcs

    # remove sources
    dataset.image.data = np.zeros_like(dataset.image.data)

    #----------------------------------------
    # define 2 spherical sky circles
    sph_circle1 = CircleSphericalSkyRegion(
        center=SkyCoord(20, 0, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))

    sph_circle2 = CircleSphericalSkyRegion(
        center=SkyCoord(50, 45, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))

    # define skycoords
    lon = np.arange(-180, 181, 10)
    lat = np.arange(-90, 91, 10)
    coords = np.array(np.meshgrid(lon, lat)).T.reshape(-1, 2)
    skycoords = SkyCoord(coords, unit='deg', frame='galactic')

    #----------------------------------------
    # get events in AND and XOR
    sph_compound_and = sph_circle1 & sph_circle2
    sph_compound_xor = sph_circle1 ^ sph_circle2

    sph_mask_and = sph_compound_and.contains(skycoords)
    sph_skycoords_and = skycoords[sph_mask_and]
    sph_mask_xor = sph_compound_xor.contains(skycoords)
    sph_skycoords_xor = skycoords[sph_mask_xor]

    # plot
    fig = plt.figure()
    fig.set_size_inches(7,3.5)
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs, aspect='equal')

    ax.scatter(skycoords.l.value, skycoords.b.value, label='all',
            transform=ax.get_transform('galactic'))
    ax.scatter(sph_skycoords_xor.l.value, sph_skycoords_xor.b.value, color='orange',
            label='xor', transform=ax.get_transform('galactic'))
    ax.scatter(sph_skycoords_and.l.value, sph_skycoords_and.b.value, color='magenta',
            label='and', transform=ax.get_transform('galactic'))

    boundary_kwargs = dict(
        include_boundary_distortions=True, discretize_kwargs={"n_points":1000}
    )
    sph_circle1.to_pixel(wcs=wcs,**boundary_kwargs).plot(ax=ax, edgecolor='green', facecolor='none',
                                alpha=0.8, lw=3)
    sph_circle2.to_pixel(wcs=wcs,**boundary_kwargs).plot(ax=ax, edgecolor='red', facecolor='none',
                                alpha=0.8, lw=3)

    ax.legend(loc='lower right')

    ax.set_xlim(-0.5, dataset.config['shape'][1] - 0.5)
    ax.set_ylim(-0.5, dataset.config['shape'][0] - 0.5)
    ax.set_title("Spherical SkyRegions")
