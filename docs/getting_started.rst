.. include:: references.txt

.. _gs:

***************
Getting started
***************

.. warning::
    This ``regions`` package is in a very early stage of development.
    It is not feature complete or API stable!
    That said, please have a look and try to use it for your applications.
    Feedback and contributions welcome!

.. _gs-intro:

Introduction
============

The `regions` package provides classes representing sky regions,
for example `~regions.CircleSkyRegion`, as well as "pixel regions",
for example `~regions.CirclePixelRegion`.

Sky regions are independent of an image. Pixel regions are connected to an image.
To transform between sky and pixel regions, a word coordinate system
(represented by a `~astropy.wcs.WCS` object) is used.

Some functions for region-based calculations (e.g. filtering a table of sky or
pixel positions) as well as functions for region serialisation (e.g. to and from
ds9 region string format) are available.

.. _gs-ds:

Dataset
=======

Throughout this tutorial, we will be working with the same example dataset and assume that
you have run `~regions.make_example_dataset` to create example ``dataset`` and ``wcs`` objects like this:

.. code-block:: python

    >>> from regions import make_example_dataset
    >>> dataset = make_example_dataset(data='simulated')
    >>> wcs = dataset.wcs

For image examples, we will use the ``wcs`` and ``image`` attributes.
For example positions, we will use the ``source_table`` and ``event_table`` attributes.

In your own analyses, you will usually load image data and WCS objects from file or
compute them with a Python script. We don't do this here, because we wanted to make
this tutorial independent of any example data files, to help you get started quickly.

Also, this example image was created to illustrate some of the key issues when working with
regions. It represents an all-sky image in Aitoff (``AIT``) projection, which means that
there are pixels at the edge of the iamge that don't correspond to positions on the sky.
And the pixels are huge, roughly 10 deg time 10 deg, which means that there are well-visible
differences between sky and pixel regions, caused by the WCS projection.

Before we start diving into coding with regions, here's an image that illustrates our
example counts image, with source positions and a few regions overplotted:

.. plot:: plot_example.py
    :include-source: false

.. _gs-coord:

Coordinates
===========

This regions package uses `astropy.coordinates.SkyCoord` objects to represent sky
coordinates.

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

`~astropy.coordinates.SkyCoord` is a very powerful and complex class (different representations,
a coordinate transformation tree) that is documented extensively in the `astropy.coordinates` docs.

In contrast, `~regions.PixCoord` is a small and simple helper class. The pixel coordinate is always
represented as cartesian coordinate data members ``x`` and ``y``. A pixel coordinate doesn't have a
frame and is not connected to the `astropy.coordinates` transformation tree.

For a given image, represented by a `~astropy.wcs.WCS` object, it's easy to transform back and
forth between sky and pixel coordinates:

.. code-block:: python

    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> pixcoord = PixCoord.from_sky(skycoord=skycoord, wcs=wcs)
    >>> pixcoord
    PixCoord(x=146.2575703393558, y=131.5998051082584)
    >>> pixcoord.to_sky(wcs=wcs)
    <SkyCoord (Galactic): (l, b) in deg
        (42.0, 43.0)>

This is an object-oriented thin wrapper around the functionality provided by `~astropy.wcs.WCS`
and `astropy.wcs.utils`.

It is possible to create `~astropy.coordinates.SkyCoord` and `~regions.PixCoord` objects
that represent arrays of pixel coordinates, and operations like transforming between sky-
and pixel or region containment checks work as expected (i.e. return arrays of the same
shape as the inputs, and perform operations on array entries independently.

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


To represent angles both on the sky and in an image, `~astropy.coordinates.Angle` objects
or `~astropy.units.Quantity` objects with angular units can be used.

.. _gs-sky:

Sky regions
===========

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

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list like this:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*SkyRegion?

.. _gs-pix:

Pixel regions
=============

In the previous section we introduced sky regions, which represent (surprise) regions on the sky sphere,
and are independent of an image or projection of the sky.

Sometimes you need regions connected to a sky image, where coordinates are given in cartesian pixel coordinates.
For those, there's a `~regions.PixCoord` class to represent a point, and a set of "pixel region" classes.
One example is `~regions.CirclePixelRegion`:

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

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list like this:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*PixelRegion?

.. _gs-shapes:

Shapes
======

This section shows one example how to construct a region for each shape that's currently supported.

* `~regions.CircleSkyRegion` and `~regions.CirclePixelRegion`

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import PixCoord, CircleSkyRegion, CirclePixelRegion

    circle_sky = CircleSkyRegion(
        center=SkyCoord(42, 43, unit='deg'),
        radius=Angle(3, 'deg'),
    )
    circle_pix = CirclePixelRegion(
        center=PixCoord(x=42, y=43),
        radius=4.2,
    )

* `~regions.EllipseSkyRegion` and `~regions.EllipsePixelRegion`

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import PixCoord, EllipseSkyRegion, EllipsePixelRegion

    ellipse_sky = EllipseSkyRegion(
        center=SkyCoord(42, 43, unit='deg'),
        minor=Angle(3, 'deg'),
        major=Angle(4, 'deg'),
        angle=Angle(5, 'deg'),
    )
    ellipse_pix = EllipsePixelRegion(
        center=PixCoord(x=42, y=43),
        minor=3,
        major=4,
        angle=Angle(5, 'deg'),
    )

* `~regions.PointSkyRegion` and `~regions.PointPixelRegion`

.. code-block:: python

    from astropy.coordinates import SkyCoord
    from regions import PixCoord, PointSkyRegion, PointPixelRegion

    point_sky = PointSkyRegion(
        center=SkyCoord(42, 43, unit='deg'),
    )
    point_pix = PointPixelRegion(
        center=PixCoord(x=42, y=43),
    )

* `~regions.PolygonSkyRegion` and `~regions.PolygonPixelRegion`

.. code-block:: python

    from astropy.coordinates import SkyCoord
    from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion

    polygon_sky = PolygonSkyRegion(
        vertices=SkyCoord([1, 2, 2], [1, 1, 2], unit='deg'),
    )
    polygon_pix = PolygonPixelRegion(
        vertices=PixCoord(x=[1, 2, 2], y=[1, 1, 2]),
    )

* `~regions.RectangleSkyRegion` and `~regions.RectanglePixelRegion`

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    from regions import PixCoord, RectangleSkyRegion, RectanglePixelRegion

    rectangle_sky = RectangleSkyRegion(
        center=SkyCoord(42, 43, unit='deg'),
        width=Angle(3, 'deg'),
        height=Angle(4, 'deg'),
        angle=Angle(5, 'deg'),
    )
    rectangle_pix = RectanglePixelRegion(
        center=PixCoord(x=42, y=43),
        width=3,
        height=4,
        angle=Angle(5, 'deg'),
    )

.. _gs-poly:

Polygons
========

Polygons are the most versatile region. Any region can be approximated as a polygon.

TODO: explain how polygons are implemented and special polygon methods, e.g. how to
obtain a polygon approximation for any shape.
This is not available yet, for now see `spherical_geometry`_ for spherical polygons
and `shapely`_ for pixel polygons.


.. _gs-wcs:

Transformations
===============

In the last two sections, we talked about how for every region shape (e.g. circle),
there's two classes, one representing "sky regions" and another representing "pixel regions"
on a given image.

A key feature of the regions package is that, for a given image, more precisely a given
`~astropy.wcs.WCS` object, it is possible to convert back and forth between sky and image
regions.


With this `wcs` object, it's possible to transform back and forth between sky and pixel regions.
As an example, let's use this sky circle:

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

.. _gs-contain:

Containment
===========

Let's continue with the ``sky_reg`` and ``pix_reg`` objects defined in the previous section:

.. code-block:: python

    >>> print(sky_reg)
    CircleSkyRegion
    center:<SkyCoord (ICRS): (ra, dec) in deg
        (50.0, 10.0)>
    radius:30.0 deg
    >>> print(pix_reg)
    CirclePixelRegion
    center: PixCoord(x=55.35205711214607, y=40.0958313892697)
    radius: 36.93290808340659


To test if a given point is inside or outside the regions, the Python ``in`` operator
can be called, which calls the special ``__contains__`` method defined on the region classes:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord
    >>> SkyCoord(50, 10, unit='deg') in sky_reg
    True
    >>> SkyCoord(50, 60, unit='deg') in sky_reg
    False
    >>> PixCoord(29, 3) in pix_reg
    True
    >>> PixCoord(29, 10) in pix_reg
    False

The ``in`` operator only works for scalar coordinates,
because Python requires the return value to be a scalar bool.
If you try to use ``in`` for non-scalar coordinates, you'll get a ``ValueError``:

.. code-block:: python

    >>> skycoords = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> skycoords in sky_reg
    ValueError: coord must be scalar. coord=<SkyCoord (ICRS): (ra, dec) in deg
        [(50.0, 10.0), (50.0, 60.0)]>

If you have arrays of coordinates, use the `regions.SkyRegion.contains`
or `regions.PixelRegion.contains` methods:

.. code-block:: python

    >>> skycoords = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sky_reg.contains(skycoords)
    array([ True, False], dtype=bool)
    >>> pixcoords = PixCoord.from_sky(skycoords, wcs)
    >>> pix_reg.contains(pixcoords)
    array([ True, False], dtype=bool)

.. _gs-masks:

Masks
=====

For aperture photometry, a common operation is to compute, for a given image and
region, a mask or array of pixel indices defining which pixels (in the whole
image or a minimal rectangular bounding box) are inside and outside the region.

All :class:`~regions.PixelRegion` objects have a
:meth:`~regions.PixelRegion.to_mask` method that returns a
:class:`~regions.Mask` object that contains information about whether
pixels are inside the region, and can be used to mask data arrays:

    >>> from regions.core import PixCoord
    >>> from regions.shapes.circle import CirclePixelRegion
    >>> center = PixCoord(4., 5.)
    >>> reg = CirclePixelRegion(center, 2.3411)
    >>> mask = reg.to_mask()
    >>> mask
    <regions.mask.Mask object at 0x10900b5c0>
    >>> mask.data
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  1.,  1.,  0.,  0.],
           [ 0.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 0.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 0.,  0.,  1.,  1.,  1.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]])

The mask data contains floating point that are between 0 (no overlap) and 1
(overlap). By default, this is determined by looking only at the central position
in each pixel, and::

    >>> reg.to_mask()

is equivalent to::

    >>> reg.to_mask(mode='center')

but other modes are available:

* ``mode='exact'``: the overlap is determined using the exact geometrical
  overlap between pixels and the region. This is slower than using the central
  position, but allows partial overlap to be treated correctly.

* ``mode='subpixels'``: the overlap is determined by sub-sampling the pixel
  using a grid of sub-pixels. The number of sub-pixels to use in this mode
  should be given using the ``subpixels`` argument.

Here are what the different modes look like:

.. plot::
   :include-source:

    import matplotlib.pyplot as plt
    from regions.core import PixCoord
    from regions.shapes.circle import CirclePixelRegion

    center = PixCoord(6.6, 7.2)
    reg = CirclePixelRegion(center, 5.2)

    plt.figure(figsize=(6, 6))

    mask1 = reg.to_mask(mode='center')
    plt.subplot(2, 2, 1)
    plt.title("mode='center'", size=9)
    plt.imshow(mask1.data, cmap=plt.cm.viridis,
               interpolation='nearest', origin='lower')

    mask2 = reg.to_mask(mode='exact')
    plt.subplot(2, 2, 2)
    plt.title("mode='exact'", size=9)
    plt.imshow(mask2.data, cmap=plt.cm.viridis,
               interpolation='nearest', origin='lower')

    mask3 = reg.to_mask(mode='subpixels', subpixels=3)
    plt.subplot(2, 2, 3)
    plt.title("mode='subpixels', subpixels=3", size=9)
    plt.imshow(mask3.data, cmap=plt.cm.viridis,
               interpolation='nearest', origin='lower')

    mask4 = reg.to_mask(mode='subpixels', subpixels=20)
    plt.subplot(2, 2, 4)
    plt.title("mode='subpixels', subpixels=20", size=9)
    plt.imshow(mask4.data, cmap=plt.cm.viridis,
               interpolation='nearest', origin='lower')

As we've seen above, the :class:`~regions.Mask` objects have a ``data``
attribute that contains a Numpy array with the mask values. However, if you
have for example a circular region with a radius of 3 pixels at a pixel position
of (1000, 1000), it would be inefficient to store a mask that has a size larger
than this, so instead we store the mask using the minimal array that contains
the mask, and the :class:`~regions.Mask` objects also include a ``bbox``
attribute that is a :class:`~regions.BoundingBox` object used to indicate where
the mask should be applied in an image.

:class:`~regions.Mask` objects also have a number of methods to make it
easy to use the masks with data. The :meth:`~regions.Mask.to_image` method
can be used to obtain an image of the mask in a 2D array of the given shape.
This places the mask in the correct place in the image and deals properly with
boundary effects:

.. plot::
   :include-source:

    import matplotlib.pyplot as plt
    from regions.core import PixCoord
    from regions.shapes.circle import CirclePixelRegion

    center = PixCoord(6.6, 7.2)
    reg = CirclePixelRegion(center, 5.2)

    mask = reg.to_mask(mode='exact')
    plt.figure(figsize=(4, 4))
    plt.imshow(mask.to_image((10, 10)), cmap=plt.cm.viridis,
               interpolation='nearest', origin='lower')

The :meth:`~regions.Mask.cutout` method can be used to create a cutout from
the input data over the mask bounding box, and the
:meth:`~regions.Mask.multiply` method can be used to multiply the aperture
mask with the input data to create a mask-weighted data cutout. All of these
methods properly handle the cases of partial or no overlap of the aperture mask
with the data.

These masks can be used as the building blocks for photometry, which we
demonstrate with a simple example. We start off by getting an example image:

.. plot::
   :context: reset
   :include-source:
   :align: center
   :nofigs:

    >>> from astropy.io import fits
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('photometry/M6707HH.fits')
    >>> hdu = fits.open(filename)[0]

We then define the aperture:

.. plot::
   :context:
   :include-source:
   :align: center
   :nofigs:

    >>> from regions.core import PixCoord
    >>> from regions.shapes.circle import CirclePixelRegion
    >>> center = PixCoord(158.5, 1053.5)
    >>> aperture = CirclePixelRegion(center, 4.)

We convert the aperture to a mask and extract a cutout from the data, as well
as a cutout with the data multipled by the mask:

.. plot::
   :context:
   :include-source:
   :align: center
   :nofigs:

    >>> mask = aperture.to_mask(mode='exact')
    >>> data = mask.cutout(hdu.data)
    >>> weighted_data = mask.multiply(hdu.data)

We can take a look at the results to make sure the source overlaps with the
aperture:

.. plot::
   :context:
   :include-source:
   :align: center

    >>> import matplotlib.pyplot as plt
    >>> plt.subplot(1,3,1)
    >>> plt.title("Mask", size=9)
    >>> plt.imshow(mask.data, cmap=plt.cm.viridis,
    ...            interpolation='nearest', origin='lower',
    ...            extent=mask.bbox.extent)
    >>> plt.subplot(1,3,2)
    >>> plt.title("Data cutout", size=9)
    >>> plt.imshow(data, cmap=plt.cm.viridis,
    ...            interpolation='nearest', origin='lower',
    ...            extent=mask.bbox.extent)
    >>> plt.subplot(1,3,3)
    >>> plt.title("Data cutout multiplied by mask", size=9)
    >>> plt.imshow(weighted_data, cmap=plt.cm.viridis,
    ...            interpolation='nearest', origin='lower',
    ...            extent=mask.bbox.extent)

We can also use the ``Mask.bbox`` attribute to look at the extent
of the mask in the image:

.. plot::
   :context:
   :include-source:
   :align: center

    >>> ax = plt.subplot(1, 1, 1)
    >>> ax.imshow(hdu.data, cmap=plt.cm.viridis,
    ...            interpolation='nearest', origin='lower')
    >>> ax.add_patch(mask.bbox.as_patch(facecolor='none', edgecolor='white'))
    >>> ax.add_patch(aperture.as_patch(facecolor='none', edgecolor='orange'))
    >>> ax.set_xlim(120, 180)
    >>> ax.set_ylim(1020, 1080)

Finally, we can use the mask and data values to compute weighted statistics:

.. plot::
   :context:
   :include-source:
   :align: center
   :nofigs:

    >>> np.average(data, weights=mask)
    9364.0126748880211

.. _gs-compound:

Compounds
=========

There's a few ways to combine any two `~regions.Region` objects into a compound region,
i.e. a `~regions.CompoundPixelRegion` or `~regions.CompoundSkyRegion` object.

* The ``&`` operator calls the ``__and__`` method which calls the :meth:`~regions.Region.intersection` method
  to create an intersection compound region.
* The ``|`` operator calls the ``__or__`` method which calls the :meth:`~regions.Region.union` method
  to create a union compound region.
* The ``^`` operator calls the ``__xor__`` method which calls the :meth:`~regions.Region.symmetric_difference` method
  to create a symmetric difference compound region.

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion
    >>> circle1 = CircleSkyRegion(
    ...    center=SkyCoord(1,2, unit='deg', frame='galactic'),
    ...    radius=Angle('5 deg')
    ...    )
    >>> circle2 = CircleSkyRegion(
    ...    center=SkyCoord(-4,3, unit='deg', frame='galactic'),
    ...    radius=Angle('3 deg'),
    ...    )
    >>> type(circle1 & circle2)
    regions.core.compound.CompoundSkyRegion
    >>> print(circle1 ^ circle2)
    (CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
    (1.0, 2.0)>
    radius: 5.0 deg
    <built-in function xor>
    CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
    (356.0, 3.0)>
    radius: 3.0 deg)

.. plot:: plot_compound.py
    :include-source: false

.. _gs-lists:

Lists
=====

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

.. _gs-shapely:

Shapely
=======

The `shapely`_ Python package is a generic package for the manipulation and analysis of geometric objects in the
Cartesian plane. Concerning regions in the cartesian plane, it is more feature-complete, powerful and optimized
than this ``regions`` package. It doesn't do everything astronomers need though, e.g. no sky regions,
no use of Astropy classes like `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.Angle`, and no
region serialisation with the formats astronomers use (such as e.g. ds9 regions).

`~regions.PixelRegion` classes provide a method :meth:`~regions.PixelRegion.to_shapely` that allows creation
of Shapely shape objects. At the moment there is no ``from_shapely`` method to convert Shapely objects
back to ``regions`` objects. The future of the use of Shapely in ``regions`` is currently unclear, some options are:

1. Add ``from_shapely`` and use it to implement e.g. `~regions.PolygonPixelRegion` operations
   can "discretization" of other shapes to polygons.
   This would make Shapely a required dependency to work with polygons.
2. Keep ``to_shapely`` for the (small?) fraction of users that want to do this,
   but don't expand or use it inside the ``regions`` package  to avoid the extra heavy dependency.
3. Remove the use of Shapely completely from the API unless good use cases demonstrating a need come up.

Here's an example how to create a Shapely object and do something that's not implemented in ``regions``,
namely to buffer a rectangle, resulting in a polygon.

.. code-block:: python

    # TODO: RectanglePixelRegion isn't implemented yet, this doesn't work yet.
    from regions import RectanglePixelRegion
    region = RectanglePixelRegion(center=(3, 2), width=2, height=1)
    shape = region.to_shapely()
    # `shape` is a `shapely.geometry.polygon.Polygon` object
    shape2 = shape.buffer(distance=3)
    # `shape2` is a polygon that's buffered by 3 pixels compared to `shape`

.. _gs-ds9:

DS9
===

The regions package provides functions to serialise and de-serialise Python lists of
`~regions.Region` objects to DS9 region strings: `~regions.ds9_objects_to_string`
and `~regions.ds9_string_to_objects`.

.. code-block:: python

    >>> from regions import ds9_objects_to_string, ds9_string_to_objects
    >>> reg_string = 'galactic\ncircle(42,43,3)'
    >>> regions = ds9_string_to_objects(reg_string)
    >>> regions[0]
    CircleSkyRegion
    center: <Galactic Coordinate: (l, b) in deg
        (42.0, 43.0)>
    radius: 3.0 deg
    >>> ds9_objects_to_string(regions, coordsys='galactic')
    '# Region file format: DS9 astropy/regions\ngalactic\ncircle(42.0000,43.0000,3.0000)\n'

There's also `~regions.write_ds9` and `~regions.read_ds9` with write to and read from
a file in addition to doing the region serialisation and parsing.

.. code-block:: python

    >>> from regions import read_ds9, write_ds9
    >>> filename = 'ds9.reg'
    >>> write_ds9(regions, filename)
    >>> regions = read_ds9(filename)
    >>> regions
    [CircleSkyRegion
     center: <FK5 Coordinate (equinox=J2000.000): (ra, dec) in deg
         (245.3477, 24.4291)>
     radius: 3.0 deg]

Often DS9 region files contain extra information about colors or other attributes.
This information is lost when converting to `~regions.Region` objects.

To make it available, the following two functions are made available:
`~regions.ds9_string_to_region_list`, `~regions.ds9_region_list_to_objects`.
Together they make up the `~regions.ds9_string_to_objects` function, but they
expose the intermediate "region list", which contains the extra attributes.

.. code-block:: python

    >>> from regions import ds9_string_to_region_list, ds9_region_list_to_objects
    >>> reg_string = """
    ... # Region file format: DS9 version 4.1
    ... global color=green dashlist=8 3
    ... galactic
    ... circle(42,43,3) # color=pink width=3 A comment
    ... circle(99,33,1) # width=2
    ... """
    >>> region_list = ds9_string_to_region_list(reg_string)
    >>> region_list
    [('circle', [<Galactic Coordinate: (l, b) in deg
           (42.0, 43.0)>, <Quantity 3.0 deg>], {'color': 'pink',
       'include': '',
       'width': '3 A comment'}),
     ('circle', [<Galactic Coordinate: (l, b) in deg
           (99.0, 33.0)>, <Quantity 1.0 deg>], {'include': '', 'width': '2'})]
    >>> regions = ds9_region_list_to_objects(region_list)
    >>> regions
    [CircleSkyRegion
     center: <Galactic Coordinate: (l, b) in deg
         (42.0, 43.0)>
     radius: 3.0 deg, CircleSkyRegion
     center: <Galactic Coordinate: (l, b) in deg
         (99.0, 33.0)>
     radius: 1.0 deg]

.. warning::

    This is very confusing, because there are two "region lists", one with tuples
    and one with `~regions.Region` objects as input. Need to find a better API
    or at least better names.

    Also, all regions currently have ``meta`` and ``visual`` arguments for ``__init__``
    and stored as region data members. These need to be documented and tests added,
    or removed.

.. _gs-mpl:

Matplotlib plotting
===================

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

.. _gs-what-next:

What next?
==========

Congratulations, you have made it to the end of the tutorial getting started guide of the
``regions`` package.

For detailed information on some specific functionality, see the API documentation here: `regions`.
Be warned though that a lot of methods haven't been implemented yet.
If you try them, they will raise a ``NotImplementedError``.

If you have the skills and time, please head over to
https://github.com/astropy/regions
and help out with ``regions`` development.
Of course, feature requests are also welcome ... they help us prioritize development efforts.
