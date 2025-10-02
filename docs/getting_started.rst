.. include:: references.txt

.. _getting_started:

Getting Started
===============

Introduction
------------

The Regions package provides classes to represent:

* Regions defined using pixel coordinates (a "region-on-image"; e.g.,
  `~regions.CirclePixelRegion`)
* Regions defined using celestial coordinates, but still in an Euclidean
  geometry (i.e., a planar projection, as a "region-on-image";
  e.g., `~regions.CircleSkyRegion`)
* Regions defined using celestial coordinates, and with a spherical
  geometry (a "region-on-celestial-sphere"; e.g., `~regions.CircleSphericalSkyRegion`)

To transform between (planar) sky and pixel regions, a `world coordinate system
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object (e.g.,
`astropy.wcs.WCS`) is needed. To transform between spherical and planar (sky or pixel)
regions, in addition to a `wcs
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_, it is also
necessary to specify whether or not boundary distortions should be included
(capturing the projection effects inherent in projection-to-spherical
transformations, or the inverse).

Regions also provides a unified interface for reading, writing,
parsing, and serializing regions data in different formats, including
the `DS9 Region Format <http://ds9.si.edu/doc/ref/region.html>`_,
`CRTF (CASA Region Text Format)
<https://casadocs.readthedocs.io/en/stable/notebooks/image_analysis.html
#Region-File-Format>`_, and `FITS Region Binary Table
<https://fits.gsfc.nasa.gov/registry/region.html>`_ format.


.. _getting_started-coord:

Coordinates
-----------

Pixel Coordinates
~~~~~~~~~~~~~~~~~

:class:`~regions.PixCoord` objects are used to represent pixel
coordinates. Pixel coordinates are defined with a scalar or an array of
``x`` and ``y`` Cartesian coordinates:

.. code-block:: python

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
array of pixel coordinates:

.. code-block:: python

    >>> pixcoord = PixCoord(x=[0, 1], y=[2, 3])
    >>> pixcoord
    PixCoord(x=[0 1], y=[2 3])
    >>> pixcoord.x
    array([0, 1])
    >>> pixcoord.y
    array([2, 3])
    >>> pixcoord.xy
    (array([0, 1]), array([2, 3]))

Let's now create a 2D array of pixel coordinates:

.. code-block:: python

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

Let's create a single sky coordinate:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> skycoord
    <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>

Sky coordinates also support array coordinates. These work in the same
way as single-value coordinates, but they store multiple coordinates in
a single object:

.. code-block:: python

    >>> skycoord = SkyCoord(ra=[10, 11, 12], dec=[41, 42, 43], unit='deg')
    >>> skycoord
    <SkyCoord (ICRS): (ra, dec) in deg
    [(10., 41.), (11., 42.), (12., 43.)]>

To represent angles both on the sky and in an image,
`~astropy.coordinates.Angle` objects or `~astropy.units.Quantity`
objects with angular units can be used.


Pixel/Sky Coordinate Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To transform between pixel and sky coordinates, a `world coordinate system
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object (e.g.,
`astropy.wcs.WCS`) is needed.

Let's start by creating a WCS object:

.. code-block:: python

    >>> from regions import make_example_dataset
    >>> dataset = make_example_dataset(data='simulated')
    >>> wcs = dataset.wcs

Now let's use this WCS to convert between sky and pixel coordinates:

.. code-block:: python

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
regions classes, for example, `~regions.CirclePixelRegion`:

.. code-block:: python

    >>> from regions import PixCoord, CirclePixelRegion
    >>> center = PixCoord(x=42, y=43)
    >>> radius = 4.2
    >>> region = CirclePixelRegion(center, radius)

You can print the region to get some information about its properties:

.. code-block:: python

    >>> print(region)
    Region: CirclePixelRegion
    center: PixCoord(x=42, y=43)
    radius: 4.2

You can access its properties via attributes:

.. code-block:: python

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
Please note they are **not** defined as regions on the celestial sphere,
but rather are meant to represent shapes on an image ("region-on-image").
They simply use sky coordinates instead of pixel coordinates to define
their position. The remaining shape parameters are converted to pixels
using the pixel scale of the image.

Let's create a sky region:

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = Angle(3, 'deg')
    >>> region = CircleSkyRegion(center, radius)

Alternatively, one can define the radius using a
`~astropy.units.Quantity` object with angular units:

.. code-block:: python

    >>> import astropy.units as u
    >>> from regions import CircleSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = 3.0 * u.deg
    >>> region = CircleSkyRegion(center, radius)

You can print the region to get some information about its properties:

.. code-block:: python

    >>> print(region)
    Region: CircleSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg

You can access its properties via attributes:

.. code-block:: python

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

Spherical sky regions are defined using celescial coordinates,
and **are** defined as regions on the celestial sphere
("regions-on-celestial-sphere", in contrast to the planar Sky Regions).
In order to transform between spherical and planar ("region-on-image") regions,
the planar projection (encoded in a `world coordinate system
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object; e.g.,
`astropy.wcs.WCS`) must be specified, along with a specification of
whether or not boundary distortions should be included.
These distortions (implemented through discrete boundary sampling)
capture the impact of the spherical-to-planar (or vice versa) projection.
However, it is possible to ignore these distortions (e.g.,
transforming a spherical circle to a planar circle).

Spherical sky regions are created using celestial coordinates (as
`~astropy.coordinates.SkyCoord`) and angular distances,
for instance specified as

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSphericalSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = Angle(3, 'deg')
    >>> region = CircleSphericalSkyRegion(center, radius)

Alternatively, one can define the radius using a
`~astropy.units.Quantity` object with angular units:

.. code-block:: python

    >>> import astropy.units as u
    >>> from regions import CircleSphericalSkyRegion
    >>> center = SkyCoord(42, 43, unit='deg')
    >>> radius = 3.0 * u.deg
    >>> region = CircleSphericalSkyRegion(center, radius)

You can print the region to get some information about its properties:

.. code-block:: python

    >>> print(region)
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg

You can access its properties via attributes:

.. code-block:: python

   >>> region.center
    <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
   >>> region.radius
   <Quantity 3. deg>

See the :ref:`shapes` documentation for the complete list of pixel-based
regions and to learn more about :class:`~regions.Region` objects and
their capabilities.
