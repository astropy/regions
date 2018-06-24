.. include:: references.txt

.. _sh:

.. _sh-shapes:

Shapes
======

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

.. .. _sh-poly:
..
.. Polygons
.. --------
..
.. Polygons are the most versatile region, since any region can be approximated
.. as a polygon.
..
.. TODO: explain how polygons are implemented and special polygon methods,
.. e.g. how to obtain a polygon approximation for any shape.
.. This is not available yet, for now see `spherical_geometry`_
.. for spherical polygons and `Shapely`_ for pixel polygons.

.. _sh-wcs:

Transformations
---------------

In the last two sections, we talked about how for every region shape (e.g.
circle), there are two classes, one representing "sky regions" and another
representing "pixel regions" on a given image.

A key feature of the regions package is that, for a given image, more precisely
a given `~astropy.wcs.WCS` object, it is possible to convert back and forth
between sky and image regions.

As an example, let's use this :class:`~regions.CircleSkyRegion`, a sky circle region:

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion

    >>> center = SkyCoord(50, 10, unit='deg')
    >>> radius = Angle(30, 'deg')
    >>> sky_reg = CircleSkyRegion(center, radius)

To convert it to a :class:`~regions.PixelRegion`, call the
:meth:`~regions.SkyRegion.to_pixel` method:

.. code-block:: python

    >>> pix_reg = sky_reg.to_pixel(wcs)
    >>> print(pix_reg)
    CirclePixelRegion
    center: PixCoord(x=55.35205711214607, y=40.0958313892697)
    radius: 36.93290808340659

Also to convert a :class:`~regions.PixelRegion` to a
:class:`~regions.SkyRegion`, call the :meth:`~regions.PixelRegion.to_sky` method:

.. code-block:: python

    >>> sky_reg = pix_reg.to_sky(wcs)
    >>> print(sky_reg)
    Region: CircleSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (50., 10.)>
    radius: 30.0 deg

.. _sh-meta:

Meta Data
---------

A :class:`~regions.Region` has ``meta`` and ``visual`` attribute stores the meta
data of the region. Since, this package supports various file formats
it is necessary to handle the meta attributes supported by them. To handle them
efficiently there are :class:`~regions.RegionMeta` and
:class:`~regions.RegionVisual` for meta and visual attributes.

The meta attribute provides additional data about regions such labels, tags,
comments, name, options io(DS9/CRTF) specific.
These classes , for now, just checks whether the key is valid or not.

The valid keys for RegionMeta class are:

1. ``label``:
2. ``symbol``/``point``: CRTF, DS9 (Symbol for which a point region is described)

3. ``include``: CRTF, DS9 (Region inclusion)
    - Possible Value: True, False
    - Ex: meta['include'] = True
4. ``frame``: CRTF (Frequency/Velocity Axis)
    - Possible values: 'REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB'
    - Default: image value
    - Ex: meta['frame'] = 'TOPO'

5. ``range``: CRTF (Frequency/Velocity Range)
    - Possible units: GHz, MHz, kHz, km/s, Hz, channel, chan (=channel)
    - Default: image range
    - Format: [min, max]
    - Ex: meta['range'] = [-320 * u.m/u.sec, -330 * u.m/u.s]

6. ``veltype``: CRTF (Velocity Calculation)
    - Possible values: 'RADIO', 'OPTICAL', 'Z', 'BETA', 'GAMMA'
    - Default: image value
    - Ex: meta['veltype'] = 'RADIO'

7. ``restfreq``: CRTF (Rest Frequency)
    - Possible values: `~astropy.units.Quantity` object
    - Default: image value
    - Ex: meta['restfreq'] = Quantity("1.42GHz")

8. ``corr``: CRTF (Correlational Axis)
    - Possible values: 'I', 'Q', 'U', 'V', 'RR', 'RL', 'LR', 'LL', 'XX', 'XY',
    'YX', 'YY', 'RX', 'RY', 'LX', 'LY', 'XR', 'XL', 'YR', 'YL', 'PP', 'PQ', 'QP'
    , 'QQ', 'RCircular', 'LCircular', 'Linear', 'Ptotal', 'Plinear', 'PFtotal',
     'PFlinear', 'Pangle
    - Default: all planes present in image
    - corr=[X, Y]

9. ``comment``:
10. ``coord``:
11. ``line``:
12. ``name``:
13. ``select``:
14. ``highlite``:
15. ``fixed``:
16. ``edit``:
17. ``move``:
18. ``rotate``:
19. ``delete``:
20. ``source``:
21. ``background``:
``tag``:
The visual attributes are meta data meant to be used to visualize regions, especially
used by plotting libraries such as `Matplotlib`_ in future.

The valid keys for RegionVisual class are:

1. ``color``:
2. ``dash``:
3. ``font``:
4. ``dashlist``:
5. ``symsize``:
6. ``symthick``:
7. ``fontsize``:
8. ``fontstyle``:
9. ``usetex``:
10. ``labelpos``:
11. ``labeloff``:
12. ``linewidth``:
13. ``linestyle``:
14. ``fill``:
15. ``line``:

.. _sh-lists:

Lists
-----

A `~regions.Region` object can only represent one region,
an array (a.k.a. vector or list) of regions.

This is in contrast to the aperture classes in `photutils` like
:class:`~photutils.CircularAperture` that do allow the ``positions``
(but usually not the other parameters) to be arrays:

.. doctest-skip::

    >>> from photutils import CircularAperture
    >>> positions = [(1, 2), (3, 4)]
    >>> apertures = CircularAperture(positions, r=4.2)

To represent lists of `~regions.Region` objects, you can store them in
Python lists (or other containers, but lists are the most common).
To create many similar regions or process many regions you can use for loops or
list comprehensions.

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
