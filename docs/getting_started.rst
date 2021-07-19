.. include:: references.txt

.. _getting_started:

Getting Started
===============

Introduction
------------

The Regions package provides classes to represent:

* Regions defined using pixel coordinates (e.g.,
  `~regions.CirclePixelRegion`)
* Regions defined using celestial coordinates, but still in an Euclidean
  geometry (e.g., `~regions.CircleSkyRegion`)

To transform between sky and pixel regions, a `world coordinate system
<https://docs.astropy.org/en/stable/wcs/wcsapi.html>`_ object (e.g.,
`astropy.wcs.WCS`) is needed.

Regions also provides a unified interface for reading,
writing, parsing, and serializing regions data in different
formats, including `CRTF (CASA Region Text Format)
<https://casa.nrao.edu/casadocs/casa-6.1.0/imaging/image-analysis/region
-file-format>`_, `DS9 Region Format
<http://ds9.si.edu/doc/ref/region.html>`_, and `FITS Region Binary Table
<https://fits.gsfc.nasa.gov/registry/region.html>`_.

Example Dataset
---------------

Throughout the documentation, we will be working with the same example
dataset and assume that you have run `~regions.make_example_dataset` to
create example ``dataset`` and ``wcs`` objects like this:

.. code-block:: python

    >>> from regions import make_example_dataset
    >>> dataset = make_example_dataset(data='simulated')
    >>> wcs = dataset.wcs

For image examples, we will use the ``wcs`` and ``image`` attributes.
For example positions, we will use the ``source_table`` and
``event_table`` attributes.

In your own analyses, you will usually load image data and WCS objects
from file or compute them with a Python script. We don't do this here,
because we wanted to make this tutorial independent of any example data
files, to help you get started quickly.

The example image was created to illustrate some of the key issues
when working with regions. It represents an all-sky image in Aitoff
(``AIT``) projection, which means that there are pixels at the edge of
the image that don't correspond to positions on the sky. And the pixels
are huge, roughly :math:`10 \times 10 \deg`, which means that there are
well-visible differences between sky and pixel regions, caused by the
WCS projection.

Before we start diving into coding with regions, here's an image that
illustrates our example counts image, with source positions and a few
regions overplotted:

.. plot:: plot_example.py
    :include-source: false


.. _getting_started-coord:

Coordinates
-----------

This regions package uses :class:`~astropy.coordinates.SkyCoord` objects
to represent sky coordinates.

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> skycoord
    <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>

To represent pixel coordinates, :class:`~regions.PixCoord` objects are
used.

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

`~astropy.coordinates.SkyCoord` is a very powerful and complex class
(different representations, a coordinate transformation tree) that is
documented extensively in the `astropy.coordinates` docs.

In contrast, `~regions.PixCoord` is a small and simple helper class.
The pixel coordinate is always represented as cartesian coordinate data
members ``x`` and ``y``. A pixel coordinate doesn't have a frame and is
not connected to the `astropy.coordinates` transformation tree.

A WCS object is used to transform back and forth between sky and pixel
coordinates:

.. code-block:: python

    >>> skycoord = SkyCoord(42, 43, unit='deg', frame='galactic')
    >>> pixcoord = PixCoord.from_sky(skycoord=skycoord, wcs=wcs)
    >>> pixcoord   # doctest: +FLOAT_CMP
    PixCoord(x=146.2575703393558, y=131.5998051082584)
    >>> pixcoord.to_sky(wcs=wcs)
    <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>

It is also possible to create `~astropy.coordinates.SkyCoord` and
`~regions.PixCoord` objects that represent arrays of pixel coordinates.
With such coordinates, operations like transforming between sky and
pixel or region containment checks work as expected (i.e., they return
arrays of the same shape as the inputs, and perform operations on array
entries independently).

.. code-block:: python

    # One-dimensional array of pixel coordinates
    >>> pixcoord = PixCoord(x=[0, 1], y=[2, 3])
    >>> pixcoord
    PixCoord(x=[0 1], y=[2 3])

    # Two-dimensional array pixel coordinates:
    >>> pixcoord = PixCoord(x=[[1, 2, 3], [4, 5, 6]],
    ...                     y=[[11, 12, 13], [14, 15, 16]])
    >>> print(pixcoord)
    PixCoord(x=[[1 2 3]
     [4 5 6]], y=[[11 12 13]
     [14 15 16]])

To represent angles both on the sky and in an image,
`~astropy.coordinates.Angle` objects or `~astropy.units.Quantity`
objects with angular units can be used.


Sky regions
-----------

Sky regions are regions that are defined using celestial coordinates.
Note that these are **not** defined as regions on the celestial sphere,
but rather are meant to represent shapes on an image, but simply defined
using celestial coordinates as opposed to pixel coordinates.

This is how to create a sky region:

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

You can print the regions to get some info about its properties:

.. code-block:: python

    >>> print(region)
    Region: CircleSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list using:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*SkyRegion?


Pixel regions
-------------

In some cases you might instead want to directly represent a region in
pixel coordinates. For those, there's a `~regions.PixCoord` class to
represent an (x, y) pixel coordinate and a set of "pixel-based region"
classes. One example is `~regions.CirclePixelRegion`:

.. code-block:: python

    >>> from regions import PixCoord, CirclePixelRegion

    >>> center = PixCoord(x=42, y=43)
    >>> radius = 4.2
    >>> region = CirclePixelRegion(center, radius)

You can print the regions to get some info about its properties:

.. code-block:: python

    >>> print(region)
    Region: CirclePixelRegion
    center: PixCoord(x=42, y=43)
    radius: 4.2

To see a list of all available sky regions, you can go to the API docs
or in IPython print the list using:

.. code-block:: none

    In [1]: import regions
    In [2]: regions.*PixelRegion?

To learn more about :class:`~regions.Region` objects and their
capabilities see the :ref:`shapes` documentation.
