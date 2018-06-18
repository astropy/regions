.. include:: references.txt

.. _gs:

Getting started
===============

.. _gs-intro:

Introduction
------------

The `regions` package provides (or plans to provide) classes to represent:

* Regions defined using pixel coordinates (`~regions.CirclePixelRegion`)
* Regions defined using celestial coordinates, but still in an Euclidean geometry (`~regions.CircleSkyRegion`)
* Regions defined on the celestial sphere (no examples yet)

To transform between sky and pixel regions, a world coordinate system
(represented by a `~astropy.wcs.WCS` object) is needed. Some functions for
region-based calculations (e.g. filtering a table of sky or pixel positions) as
well as functions for region serialisation (e.g. to and from ds9/crtf region string
format) are available.

.. _gs-ds:

Dataset
-------

Throughout the documentation, we will be working with the same example dataset
and assume that you have run `~regions.make_example_dataset` to create example
``dataset`` and ``wcs`` objects like this:

.. code-block:: python

    >>> from regions import make_example_dataset
    >>> dataset = make_example_dataset(data='simulated')
    >>> wcs = dataset.wcs

For image examples, we will use the ``wcs`` and ``image`` attributes. For
example positions, we will use the ``source_table`` and ``event_table``
attributes.

In your own analyses, you will usually load image data and WCS objects from file
or compute them with a Python script. We don't do this here, because we wanted
to make this tutorial independent of any example data files, to help you get
started quickly.

Also, this example image was created to illustrate some of the key issues when
working with regions. It represents an all-sky image in Aitoff (``AIT``)
projection, which means that there are pixels at the edge of the image that
don't correspond to positions on the sky. And the pixels are huge, roughly 10
deg time 10 deg, which means that there are well-visible differences between sky
and pixel regions, caused by the WCS projection.

Before we start diving into coding with regions, here's an image that
illustrates our example counts image, with source positions and a few regions
overplotted:

.. plot:: plot_example.py
    :include-source: false

.. _gs-coord:

Coordinates
-----------

This regions package uses :class:`~astropy.coordinates.SkyCoord` objects to
represent sky coordinates.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> skycoord
    <SkyCoord (Galactic): (l, b) in deg
    (42.0, 43.0)>

To represent pixel coordinates, `regions.PixCoord` objects are used.

.. code-block:: python

    >>> from regions import PixCoord
    >>> pixcoord = PixCoord(x=42, y=43)
    >>> pixcoord
    PixCoord(x=42, y=43)
    >>> pixcoord.x
    42
    >>> pixcoord.y
    43

`~astropy.coordinates.SkyCoord` is a very powerful and complex class (different
representations, a coordinate transformation tree) that is documented
extensively in the `astropy.coordinates` docs.

In contrast, `~regions.PixCoord` is a small and simple helper class. The pixel
coordinate is always represented as cartesian coordinate data members ``x`` and
``y``. A pixel coordinate doesn't have a frame and is not connected to the
`astropy.coordinates` transformation tree.

For a given image, represented by a `~astropy.wcs.WCS` object, it's easy to
transform back and forth between sky and pixel coordinates:

.. code-block:: python

    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> pixcoord = PixCoord.from_sky(skycoord=skycoord, wcs=wcs)
    >>> pixcoord
    PixCoord(x=146.2575703393558, y=131.5998051082584)
    >>> pixcoord.to_sky(wcs=wcs)
    <SkyCoord (Galactic): (l, b) in deg
        (42.0, 43.0)>

This is an object-oriented thin wrapper around the functionality provided by
`~astropy.wcs.WCS` and `astropy.wcs.utils`.

It is possible to create `~astropy.coordinates.SkyCoord` and `~regions.PixCoord`
objects that represent arrays of pixel coordinates, and operations like
transforming between sky- and pixel or region containment checks work as
expected (i.e. return arrays of the same shape as the inputs, and perform
operations on array entries independently.

.. code-block:: python

    # One-dimensional array of pixel coordinates
    >>> pixcoord = PixCoord(x=[0, 1], y=[2, 3])
    >>> pixcoord
    PixCoord(x=[0 1], y=[2 3])

    # Two-dimensional array pixel coordinates:
    >>> pixcoord = PixCoord(
    ...     x=[[1, 2, 3], [4, 5, 6]],
    ...     y=[[11, 12, 13], [14, 15, 16]]
    ... )
    >>> print(pixcoord)
    PixCoord(x=[[1 2 3]
     [4 5 6]], y=[[11 12 13]
     [14 15 16]])


To represent angles both on the sky and in an image,
`~astropy.coordinates.Angle` objects or `~astropy.units.Quantity` objects with
angular units can be used.

.. _gs-sky:

Sky regions
-----------

Sky regions are regions that are defined using celestial coordinates. Note that
these are **not** defined as regions on the celestial sphere, but rather are
meant to represent shapes on an image, but simply defined using celestial
coordinates as opposed to pixel coordinates.

This is how to create a sky region:

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import CircleSkyRegion

    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)

You can print the regions to get some info about its properties:

.. code-block:: python

    >>> print(region)
    CircleSkyRegion
    center:<SkyCoord (ICRS): (ra, dec) in deg
        (42.0, 43.0)>
    radius:3.0 deg

To see a list of all available sky regions, you can go to the API docs or in
IPython print the list using:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*SkyRegion?

.. _gs-pix:

Pixel regions
-------------

In some cases you might instead want to directly represent a region in pixel
coordinates. For those, there's a `~regions.PixCoord` class to represent a
point, and a set of "pixel region" classes. One example is
`~regions.CirclePixelRegion`:

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import PixCoord, CirclePixelRegion

    center = PixCoord(x=42, y=43)
    radius = 4.2
    region = CirclePixelRegion(center, radius)

You can print the regions to get some info about its properties:

.. code-block:: python

    >>> print(region)
    CirclePixelRegion
    center: PixCoord(x=42, y=43)
    radius: 4.2

To see a list of all available sky regions, you can go to the API docs or in
IPython print the list using:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*PixelRegion?

.. _gs-shapes:

Shapes
------

This section shows one example how to construct a region for each shape that's
currently supported.

* `~regions.CircleSkyRegion` and `~regions.CirclePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import PixCoord, CircleSkyRegion, CirclePixelRegion

    >>> circle_sky = CircleSkyRegion(center=SkyCoord(42, 43, unit='deg'),
    ...                              radius=Angle(3, 'deg'))
    >>> circle_pix = CirclePixelRegion(center=PixCoord(x=42, y=43),
    ...                                radius=4.2)

* `~regions.CircleAnnulusSkyRegion` and `~regions.CircleAnnulusPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import PixCoord, CircleAnnulusSkyRegion, CircleAnnulusPixelRegion

    >>> circle_annulus_sky = CircleAnnulusSkyRegion(center=SkyCoord(42, 43, unit='deg'),
    ...                               inner_radius=Angle(3, 'deg'), outer_radius=Angle(4, 'deg'))
    >>> circle_annulus_pix = CircleAnnulusPixelRegion(center=PixCoord(x=42, y=43),
    ...                                inner_radius=4.2, outer_radius=5.2)

* `~regions.EllipseSkyRegion` and `~regions.EllipsePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import PixCoord, EllipseSkyRegion, EllipsePixelRegion

    >>> ellipse_sky = EllipseSkyRegion(center=SkyCoord(42, 43, unit='deg'),
    ...                                height=Angle(3, 'deg'), width=Angle(4, 'deg'),
    ...                                angle=Angle(5, 'deg'))
    >>> ellipse_pix = EllipsePixelRegion(center=PixCoord(x=42, y=43),
    ...                                  height=4.2, width=4.2,
    ...                                  angle=Angle(5, 'deg'))

* `~regions.PointSkyRegion` and `~regions.PointPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, PointSkyRegion, PointPixelRegion

    >>> point_sky = PointSkyRegion(center=SkyCoord(42, 43, unit='deg'))
    >>> point_pix = PointPixelRegion(center=PixCoord(x=42, y=43))

* `~regions.LineSkyRegion` and `~regions.LinePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, LineSkyRegion, LinePixelRegion

    >>> line_sky = LineSkyRegion(start=SkyCoord(42, 43, unit='deg'), end=SkyCoord(52, 53, unit='deg'))
    >>> line_pix = LinePixelRegion(start=PixCoord(x=42, y=43), start=PixCoord(x=52, y=53))

* `~regions.PolygonSkyRegion` and `~regions.PolygonPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion

    >>> polygon_sky = PolygonSkyRegion(vertices=SkyCoord([1, 2, 2], [1, 1, 2], unit='deg'),
    >>> polygon_pix = PolygonPixelRegion(vertices=PixCoord(x=[1, 2, 2], y=[1, 1, 2]))

* `~regions.RectangleSkyRegion` and `~regions.RectanglePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import PixCoord, RectangleSkyRegion, RectanglePixelRegion

    >>> rectangle_sky = RectangleSkyRegion(center=SkyCoord(42, 43, unit='deg'),
    ...                                    width=Angle(3, 'deg'), height=Angle(4, 'deg'),
    ...                                    angle=Angle(5, 'deg'))
    >>> rectangle_pix = RectanglePixelRegion(center=PixCoord(x=42, y=43),
    ...                                    width=3, height=4,
    ...                                    angle=Angle(5, 'deg'))

.. .. _gs-poly:
..
.. Polygons
.. --------
..
.. Polygons are the most versatile region, since any region can be approximated as a polygon.
..
.. TODO: explain how polygons are implemented and special polygon methods, e.g. how to
.. obtain a polygon approximation for any shape.
.. This is not available yet, for now see `spherical_geometry`_ for spherical polygons
.. and `Shapely`_ for pixel polygons.

.. _gs-wcs:

Transformations
---------------

In the last two sections, we talked about how for every region shape (e.g.
circle), there are two classes, one representing "sky regions" and another
representing "pixel regions" on a given image.

A key feature of the regions package is that, for a given image, more precisely
a given `~astropy.wcs.WCS` object, it is possible to convert back and forth
between sky and image regions.

As an example, let's use this sky circle region:

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import CircleSkyRegion

    center = SkyCoord(50, 10, unit='deg')
    radius = Angle(30, 'deg')
    sky_reg = CircleSkyRegion(center, radius)

To convert it to a pixel region, call the :meth:`~regions.SkyRegion.to_pixel` method:

.. code-block:: python

    >>> pix_reg = sky_reg.to_pixel(wcs)
    >>> pix_reg
    CirclePixelRegion
    center: PixCoord(x=55.35205711214607, y=40.0958313892697)
    radius: 36.93290808340659

.. _gs-lists:

Lists
-----

A `~regions.Region` object can only represent one region, not an array (a.k.a. vector or list) of regions.

This is in contrast to the aperture classes in `photutils` like `photutils.CircularAperture` that
do allow the ``positions`` (but usually not the other parameters) to be arrays:

.. code-block:: python

    from photutils import CircularAperture
    positions = [(1, 2), (3, 4)]
    apertures = CircularAperture(positions, r=4.2)

To represent lists of `~regions.Region` objects, you can store them in Python lists
(or other containers, but lists are the most common).
To create many similar regions or process many regions you can use for loops or list comprehensions.

.. code-block:: python

    >>> from regions import PixCoord, CirclePixelRegion
    >>> regions = [
    ...    CirclePixelRegion(center=PixCoord(x, y), radius=4.2)
    ...    for x, y in [(1, 2), (3, 4)]
    ... ]
    >>> regions
    [CirclePixelRegion
     center: PixCoord(x=1, y=2)
     radius: 4.2, CirclePixelRegion
     center: PixCoord(x=3, y=4)
     radius: 4.2]
    >>> for region in regions:
    ...    print(region.center)
    PixCoord(x=1, y=2)
    PixCoord(x=3, y=4)
    >>> [region.area for region in regions]
    [55.41769440932395, 55.41769440932395]
