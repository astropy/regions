Overview
========

Introduction
------------

The Regions package provides classes to represent:

* Regions defined using pixel coordinates (a "region-on-image"; e.g.,
  `~regions.CirclePixelRegion`)
* Regions defined using celestial coordinates, but still in an Euclidean
  geometry (i.e., a planar projection, as a "region-on-image"; e.g.,
  `~regions.CircleSkyRegion`)

* Regions defined using celestial coordinates, and with a
  spherical geometry (a "region-on-celestial-sphere"; e.g.,
  `~regions.CircleSphericalSkyRegion`)

To transform between (planar) sky and pixel regions, a `world coordinate
system <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_
object (e.g., FITS `~astropy.wcs.WCS`) is needed. To transform between
spherical and planar (sky or pixel) regions, in addition to a `WCS
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_, it is also
necessary to specify whether or not boundary distortions should
be included (capturing the WCS projection effects inherent in
planar-to-spherical transformations, or the inverse).

Regions also provides a unified interface for reading, writing,
parsing, and serializing regions data in different formats, including
the `DS9 Region Format <http://ds9.si.edu/doc/ref/region.html>`_,
`CRTF (CASA Region Text Format)
<https://casadocs.readthedocs.io/en/stable/notebooks/image_analysis.html
#Region-File-Format>`_, and `FITS Region Binary Table
<https://fits.gsfc.nasa.gov/registry/region.html>`_ format.

The code and issue tracker are available at the following links:

* Code: https://github.com/astropy/regions
* Issue Tracker: https://github.com/astropy/regions/issues

Like much astronomy software, Regions is an evolving package. The
developers try to maintain backwards compatibility, but at times the
API may change if there is a benefit to doing so. If there are specific
areas you think API stability is important, please let us know as part
of the development process.


.. _getting_started-coord:

Coordinates
-----------

Pixel Coordinates
~~~~~~~~~~~~~~~~~

:class:`~regions.PixCoord` objects are used to represent pixel
coordinates. Pixel coordinates are defined with a scalar or an array of
``x`` and ``y`` Cartesian coordinates::

    >>> from regions import PixCoord
    >>> pixcoord = PixCoord(x=42, y=43)
    >>> pixcoord
    PixCoord(x=42, y=43)
    >>> pixcoord.x
    42
    >>> pixcoord.y
    43
    >>> pixcoord.xy
    (42, 43)

`~regions.PixCoord` objects can also represent arrays of pixel
coordinates. These work in the same way as single-value coordinates, but
they store multiple coordinates in a single object. Let's create a 1D
array of pixel coordinates::

    >>> pixcoord = PixCoord(x=[0, 1], y=[2, 3])
    >>> pixcoord
    PixCoord(x=[0 1], y=[2 3])
    >>> pixcoord.x
    array([0, 1])
    >>> pixcoord.y
    array([2, 3])
    >>> pixcoord.xy
    (array([0, 1]), array([2, 3]))

Let's now create a 2D array of pixel coordinates::

    >>> pixcoord = PixCoord(x=[[1, 2, 3], [4, 5, 6]],
    ...                     y=[[11, 12, 13], [14, 15, 16]])
    >>> pixcoord
    PixCoord(x=[[1 2 3]
     [4 5 6]], y=[[11 12 13]
     [14 15 16]])


Sky Coordinates
~~~~~~~~~~~~~~~

:class:`~astropy.coordinates.SkyCoord` objects are used to represent
sky coordinates. `~astropy.coordinates.SkyCoord` is a very powerful
class that provides a flexible interface for celestial coordinate
representation, manipulation, and transformation between systems. See
the extensive :ref:`astropy-coordinates` documentation for more details.

Let's create a single sky coordinate::

    >>> from astropy.coordinates import SkyCoord
    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> skycoord
    <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>

Sky coordinates also support array coordinates. These work in the same
way as single-value coordinates, but they store multiple coordinates in
a single object::

    >>> skycoord = SkyCoord(ra=[10, 11, 12], dec=[41, 42, 43], unit='deg')
    >>> skycoord
    <SkyCoord (ICRS): (ra, dec) in deg
    [(10., 41.), (11., 42.), (12., 43.)]>

To represent angles both on the sky and in an image,
`~astropy.coordinates.Angle` objects or `~astropy.units.Quantity`
objects with angular units can be used.


Pixel/Sky Coordinate Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To transform between pixel and sky coordinates, a `world coordinate
system <https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object
(e.g., `~astropy.wcs.WCS`) is needed.

Let's start by creating a WCS object::

    >>> from astropy.wcs import WCS
    >>> wcs = WCS(naxis=2)
    >>> wcs.wcs.crpix = (180, 90)
    >>> wcs.wcs.cdelt = (-1, 1)
    >>> wcs.wcs.crval = (0, 0)
    >>> wcs.wcs.ctype = ('GLON-AIT', 'GLAT-AIT')

Now let's use this WCS to convert between sky and pixel coordinates::

    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> pixcoord = PixCoord.from_sky(skycoord=skycoord, wcs=wcs)
    >>> pixcoord   # doctest: +FLOAT_CMP
    PixCoord(x=146.2575703393558, y=131.5998051082584)
    >>> pixcoord.to_sky(wcs=wcs)
    <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>


Pixel Regions
-------------

Pixel regions are regions that are defined using pixel coordinates and
sizes in pixels. The regions package provides a set of "pixel-based"
regions classes, for example, `~regions.CirclePixelRegion`::

    >>> from regions import PixCoord, CirclePixelRegion
    >>> center = PixCoord(x=42, y=43)
    >>> radius = 4.2
    >>> region = CirclePixelRegion(center, radius)

You can print the region to get some information about its properties::

    >>> print(region)
    Region: CirclePixelRegion
    center: PixCoord(x=42, y=43)
    radius: 4.2

You can access its properties via attributes::

   >>> region.center
   PixCoord(x=42, y=43)
   >>> region.radius
   4.2

See the :ref:`shapes` documentation for the complete list of pixel-based
regions and to learn more about :class:`~regions.Region` objects and
their capabilities.


Sky Regions
-----------

Sky regions are regions that are defined using celestial coordinates.
Please note they are **not** defined as regions on the celestial
sphere, but rather are meant to represent shapes on an image
("region-on-image"). They simply use sky coordinates instead of pixel
coordinates to define their position. The remaining shape parameters are
converted to pixels using the pixel scale of the image.

Let's create a sky region::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = Angle(3, 'deg')
    >>> region = CircleSkyRegion(center, radius)

Alternatively, one can define the radius using a
`~astropy.units.Quantity` object with angular units::

    >>> import astropy.units as u
    >>> from regions import CircleSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = 3.0 * u.deg
    >>> region = CircleSkyRegion(center, radius)

You can print the region to get some information about its properties::

    >>> print(region)
    Region: CircleSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg

You can access its properties via attributes::

   >>> region.center
    <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
   >>> region.radius
   <Quantity 3. deg>

See the :ref:`shapes` documentation for the complete list of pixel-based
regions and to learn more about :class:`~regions.Region` objects and
their capabilities.


Spherical Sky Regions
---------------------

Spherical sky regions are defined using celestial coordinates,
and **are** defined as regions on the celestial sphere
("regions-on-celestial-sphere", in contrast to the planar Sky Regions).

Spherical sky regions are created using celestial coordinates (as
`~astropy.coordinates.SkyCoord`) and angular distances, for instance
specified as::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSphericalSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = Angle(3, 'deg')
    >>> region = CircleSphericalSkyRegion(center, radius)

Alternatively, one can define the radius using a
`~astropy.units.Quantity` object with angular units::

    >>> import astropy.units as u
    >>> from regions import CircleSphericalSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = 3.0 * u.deg
    >>> region = CircleSphericalSkyRegion(center, radius)

You can print the region to get some information about its properties::

    >>> print(region)
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg

You can access its properties via attributes::

   >>> region.center
    <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
   >>> region.radius
   <Quantity 3. deg>

See the :ref:`shapes` documentation for the complete list of pixel-based
regions and to learn more about :class:`~regions.Region` objects and
their capabilities.

Spherical to planar region transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to transform between spherical and planar ("region-on-image")
regions, the planar projection (encoded in a `world coordinate system
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object; e.g.,
`~astropy.wcs.WCS`) must be specified, along with a specification
of whether or not boundary distortions should be included. These
distortions (implemented through discrete boundary sampling) capture the
impact of the spherical-to-planar (or vice versa) projection described
by the WCS. However, it is possible to ignore these distortions (e.g.,
transforming a spherical circle to a planar circle).

The example below demonstrates the difference between a spherical and
a planar circle with the same center and radius. Boundary distortions
are included when projecting the spherical circle onto the plot. In this
full-sky Aitoff projection, the distortions result in points inside the
spherical circle falling outside of the planar circle and vice versa.

.. plot::
    :include-source: false

    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.coordinates import Angle, SkyCoord
    from astropy.wcs import WCS
    from regions import CircleSkyRegion, CircleSphericalSkyRegion, PixCoord

    # Create a full-sky Aitoff WCS
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = (180, 90)
    wcs.wcs.cdelt = (-1, 1)
    wcs.wcs.crval = (0, 0)
    wcs.wcs.ctype = ('GLON-AIT', 'GLAT-AIT')
    shape = (180, 360)

    # Define skycoords, pixcoords grids
    lon = np.arange(-180, 181, 10)
    lat = np.arange(-90, 91, 10)
    coords = np.array(np.meshgrid(lon, lat)).T.reshape(-1, 2)
    skycoords = SkyCoord(coords, unit='deg', frame='galactic')
    pixcoords = PixCoord.from_sky(skycoords, wcs)

    # Define spherical & planar sky circles
    sph_circle = CircleSphericalSkyRegion(
        center=SkyCoord(50, 45, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))
    circle = CircleSkyRegion(
        center=SkyCoord(50, 45, unit='deg', frame='galactic'),
        radius=Angle('30 deg'))
    # Note: circle is equivalent to transforming from sph_circle
    # with sph_circle.to_sky(wcs=wcs, include_boundary_distortions=False)

    # Define transformed-to pixel regions
    pix_circ_distort = sph_circle.to_pixel(wcs=wcs,
                                           include_boundary_distortions=True,
                                           n_vertices=1000)
    pix_circ_nodistort = circle.to_pixel(wcs=wcs)

    # Get contained points
    distort_mask = sph_circle.contains(skycoords)
    nodistort_mask = pix_circ_nodistort.contains(pixcoords)

    both_skycoords = skycoords[distort_mask & nodistort_mask]
    distort_only_skycoords = skycoords[distort_mask & ~nodistort_mask]
    nodistort_only_skycoords = skycoords[~distort_mask & nodistort_mask]

    # Plot
    fig = plt.figure()
    fig.set_size_inches(7, 3.5)
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs, aspect='equal')

    ax.scatter(skycoords.l.value, skycoords.b.value, label='All',
               transform=ax.get_transform('galactic'), color='lightgrey')
    ax.scatter(distort_only_skycoords.l.value, distort_only_skycoords.b.value,
               color='magenta', label='Only within spherical circle',
               transform=ax.get_transform('galactic'))
    ax.scatter(nodistort_only_skycoords.l.value,
               nodistort_only_skycoords.b.value,
               color='lime', label='Only within planar circle',
               transform=ax.get_transform('galactic'))
    ax.scatter(both_skycoords.l.value, both_skycoords.b.value, color='orange',
               label='Within both', transform=ax.get_transform('galactic'))

    pix_circ_distort.plot(ax=ax, edgecolor='red', facecolor='none',
                          alpha=0.8, lw=3)

    pix_circ_nodistort.plot(ax=ax, edgecolor='green', facecolor='none',
                            alpha=0.8, lw=3)

    ax.legend(loc='lower right')

    ax.set_xlim(-0.5, shape[1] - 0.5)
    ax.set_ylim(-0.5, shape[0] - 0.5)
    ax.set_title('Spherical vs. Planar circle: '
                 'Center=(50°, 45°), radius=30°')
