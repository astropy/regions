.. include:: references.txt

.. _gs:

===============
Getting started
===============

The `regions` package is in a very early stage of development.
It is not feature complete or API stable.
Feedback and contributions welcome!
That said, here's a short tutorial showing the main features.

.. _gs-sky:

Sky regions
-----------

The `astropy.coordinates` package provides the `~astropy.coordinates.Angle`
and `~astropy.coordinates.SkyCoord` classes.

The `regions` package provides classes representing sky regions,
for example `~regions.CircleSkyRegion`.

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
    Center:<SkyCoord (ICRS): (ra, dec) in deg
        (42.0, 43.0)>
    Radius:3.0 deg

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list like this:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*SkyRegion?

.. _gs-pix:

Pixel regions
-------------

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
    <regions.shapes.circle.CirclePixelRegion object at 0x110038cc0>
    # TODO: implement `str` and update this example

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list like this:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*PixelRegion?

.. _gs-wcs:

Region transformations
----------------------

In the last two sections, we talked about how for every region shape (e.g. circle),
there's two classes, one representing "sky regions" and another representing "pixel regions"
on a given image.

A key feature of the regions package is that, for a given image, more precisely a given
`~astropy.wcs.WCS` object, it is possible to convert back and forth between sky and image
regions.

Usually you create the WCS object from the information in a FITS file.
For this tutorial, let's create an example ``WCS`` from scratch corresponding to an
image that is in Galactic coordinates and Aitoff projection (``ctype``)
has the reference point at the Galactic center on the sky (``crval = 0, 0``)
and at pixel coordinate ``crpix = 18, 9`` in the image, and has huge pixels
of roughly 10 deg x 10 deg (``cdelt = 10, 10``).

.. code-block:: python

    from astropy.wcs import WCS

    wcs = WCS(naxis=2)
    wcs.wcs.crval = 0, 0
    wcs.wcs.crpix = 18, 9
    wcs.wcs.cdelt = 10, 10
    wcs.wcs.ctype = 'GLON-AIT', 'GLAT-AIT'

    # shape = (36, 18) would give an image that covers the whole sky.


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
    >>> pix_reg.center
    array([ 29.36479429,   3.10958314])
    # TODO: center should be a PixCoord object
    >>> pix_reg.radius
    <Quantity 0.0010259141133922682 deg pix / arcsec>
    # TODO: radius should be a float or quantity with unit pix
    >>> pix_reg.radius.to('pix')
    <Quantity 3.6932908082121654 pix>

.. _gs-contain:

Containment
-----------

Let's continue with the ``sky_reg`` and ``pix_reg`` objects defined in the previous section:

.. code-block:: python

    >>> print(sky_reg)
    CircleSkyRegion
    Center:<SkyCoord (ICRS): (ra, dec) in deg
        (50.0, 10.0)>
    Radius:30.0 deg
    >>> print(pix_reg)
    # TODO: str isn't defined yet.

To test if a given point is inside or outside the regions, the Python ``in`` operator
can be called, which calls the special ``__contains__`` method defined on the region classes:

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord
    >>> SkyCoord(50, 10, unit='deg') in sky_reg
    True
    >>> SkyCoord(50, 60, unit='deg') in sky_reg
    False
    # TODO: PixCoord containment test currently doesn't work
    # AttributeError: 'PixCoord' object has no attribute 'isscalar'
    >>> PixCoord(29, 3) in pix_reg
    # Should be True
    >>> PixCoord(29, 10) in pix_reg
    # Should be False
    >>> pix_reg.contains(PixCoord(29, 3))
    # TODO: AttributeError: 'numpy.ndarray' object has no attribute 'x'

TODO: Example how to test many points, e.g. an event list or source catalog

It's not clear to me if the in operator can be used on sequence arrays, returning a bool array.
This description sounds like it should work:
https://docs.python.org/3/reference/datamodel.html#object.__contains__
Do we want to always use / recommend ``in``, or should we just use the
`regions.SkyRegion.contains` and `regions.PixelRegion.contains` methods?

Current behaviour:

.. code-block:: python

    >>> SkyCoord([50, 50], [10, 60], unit='deg') in sky_reg
    ValueError: <SkyCoord (ICRS): (ra, dec) in deg
    [(50.0, 10.0), (50.0, 60.0)]> must be scalar
    >>> sky_reg.contains(SkyCoord([50, 50], [10, 60], unit='deg'))
    array([ True, False], dtype=bool)

TODO: add pixel coordinate example


.. _gs-spatial:

Spatial filtering
-----------------

For aperture photometry, a common operation is to compute, for a given image and region,
a boolean mask or array of pixel indices defining which pixels (in the whole image or a
minimal rectangular bounding box) are inside and outside the region.

To a certain degree, such spatial filtering can be done using the methods described in the previous :ref:`gs-contain`
section. Apart from that, no high-level functionality for spatial filtering, bounding boxes or aperture photometry
is available yet.

For now, please use `photutils`_ or `pyregion`_.

(We plan to merge that functionality from ``photutils`` or ``pyregion`` into this ``regions`` package, or re-implement it.)

.. _gs-compound:

Compound regions
----------------

There's a few ways to combine any two `~regions.Region` objects into a compound region,
i.e. a `~regions.CompoundPixelRegion` or `~regions.CompoundSkyRegion` object.

* The ``&`` operator calls the ``__and__`` method which calls the :meth:`~regions.Region.intersection` method
  to create an intersection compound region.
* The ``|`` operator calls the ``__or__`` method which calls the :meth:`~regions.Region.union` method
  to create a union compound region.
* The ``^`` operator calls the ``__xor__`` method which calls the :meth:`~regions.Region.symmetric_difference` method
  to create a symmetric difference compound region.


TODO:

* A code example
* Add image illustrating the compound regions

.. _gs-shapely:

Shapely
-------

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

ds9 region strings
------------------

TODO: add examples for `regions.write_ds9`, `regions.read_ds9`, `regions.objects_to_ds9_string`
and `regions.region_list_to_objects`.

.. code-block:: python

    # Currently it doesn't work
    >>> from regions import objects_to_ds9_string
    >>> objects_to_ds9_string([sky_reg])
    # TypeError: 'float' object is not subscriptable
    >>> regions.write_ds9([sky_reg], filename='test.reg')
    # TypeError: 'float' object is not subscriptable
