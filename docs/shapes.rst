.. include:: references.txt

.. _sh:

.. _sh-shapes:

Shapes
======

This section shows one example how to construct a region for each shape that's
currently supported.

* `~regions.CircleSkyRegion` and `~regions.CirclePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, CircleSkyRegion, CirclePixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> circle_sky = CircleSkyRegion(center=center_sky, radius=3 * u.deg)
    >>> circle_pix = CirclePixelRegion(center=PixCoord(x=42, y=43),
    ...                                radius=4.2)

* `~regions.CircleAnnulusSkyRegion` and `~regions.CircleAnnulusPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, CircleAnnulusSkyRegion, CircleAnnulusPixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> circle_annulus_sky = CircleAnnulusSkyRegion(center=center_sky,
    ...                                             inner_radius=3 * u.deg,
    ...                                             outer_radius=4 * u.deg)
    >>> circle_annulus_pix = CircleAnnulusPixelRegion(center=PixCoord(x=42, y=43),
    ...                                               inner_radius=4.2,
    ...                                               outer_radius=5.2)

* `~regions.EllipseSkyRegion` and `~regions.EllipsePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, EllipseSkyRegion, EllipsePixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> ellipse_sky = EllipseSkyRegion(center=center_sky,
    ...                                height=3 * u.deg, width=3 * u.deg,
    ...                                angle=5 * u.deg)
    >>> ellipse_pix = EllipsePixelRegion(center=PixCoord(x=42, y=43),
    ...                                  height=4.2, width=4.2,
    ...                                  angle=5 * u.deg)

* `~regions.EllipseAnnulusSkyRegion` and `~regions.EllipseAnnulusPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, EllipseAnnulusSkyRegion, EllipseAnnulusPixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> ellipse_annulus_sky = EllipseAnnulusSkyRegion(center=center_sky,
    ...                                               inner_width=3 * u.deg,
    ...                                               outer_width=4 * u.deg,
    ...                                               inner_height=6 * u.deg,
    ...                                               outer_height=7 * u.deg,
    ...                                               angle=6 * u.deg)
    >>> ellipse_annulus_pix = EllipseAnnulusPixelRegion(center=PixCoord(x=42, y=43),
    ...                                                 inner_width=4.2,
    ...                                                 outer_width=5.2,
    ...                                                 inner_height=7.2,
    ...                                                 outer_height=8.2,
    ...                                                 angle=6 * u.deg)

* `~regions.PointSkyRegion` and `~regions.PointPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, PointSkyRegion, PointPixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> point_sky = PointSkyRegion(center=center_sky)
    >>> point_pix = PointPixelRegion(center=PixCoord(x=42, y=43))

* `~regions.TextSkyRegion` and `~regions.TextPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, TextSkyRegion, TextPixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> point_sky = TextSkyRegion(center=center_sky, 'Demo Text')
    >>> point_pix = TextPixelRegion(center=PixCoord(x=42, y=43), 'Demo Text')

* `~regions.LineSkyRegion` and `~regions.LinePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, LineSkyRegion, LinePixelRegion

    >>> start_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> end_sky = SkyCoord(52, 53, unit='deg', frame='fk5')
    >>> line_sky = LineSkyRegion(start=start_sky, end=end_sky)
    >>> line_pix = LinePixelRegion(start=PixCoord(x=42, y=43), end=PixCoord(x=52, y=53))

* `~regions.RectangleSkyRegion` and `~regions.RectanglePixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, RectangleSkyRegion, RectanglePixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> rectangle_sky = RectangleSkyRegion(center=center_sky,
    ...                                    width=3 * u.deg, height=4 * u.deg,
    ...                                    angle=5 * u.deg)
    >>> rectangle_pix = RectanglePixelRegion(center=PixCoord(x=42, y=43),
    ...                                      width=3, height=4,
    ...                                      angle=5 * u.deg)

* `~regions.RectangleAnnulusSkyRegion` and `~regions.RectangleAnnulusPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> from regions import PixCoord, RectangleAnnulusSkyRegion, RectangleAnnulusPixelRegion

    >>> center_sky = SkyCoord(42, 43, unit='deg', frame='fk5')
    >>> rectangle_annulus_sky = RectangleAnnulusSkyRegion(center=center_sky,
    ...                                                   inner_width=3 * u.deg,
    ...                                                   outer_width=4 * u.deg,
    ...                                                   inner_height=6 * u.deg,
    ...                                                   outer_height=7 * u.deg,
    ...                                                   angle=15 * u.deg)
    >>> rectangle_annulus_pix = RectangleAnnulusPixelRegion(center=PixCoord(x=42, y=43),
    ...                                                     inner_width=4.2,
    ...                                                     outer_width=5.2,
    ...                                                     inner_height=7.2,
    ...                                                     outer_height=8.2,
    ...                                                     angle=15 * u.deg)

* `~regions.PolygonSkyRegion` and `~regions.PolygonPixelRegion`

.. code-block:: python

    >>> from astropy.coordinates import SkyCoord
    >>> from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion

    >>> polygon_sky = PolygonSkyRegion(vertices=SkyCoord([1, 2, 2], [1, 1, 2], unit='deg', frame='fk5'))
    >>> polygon_pix = PolygonPixelRegion(vertices=PixCoord(x=[1, 2, 2], y=[1, 1, 2]))

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

A :class:`~regions.Region` has ``meta`` and ``visual`` attributes which
stores the meta data of the region. Since this package supports various file
formats it is necessary to handle the meta attributes supported by them.
To handle them there are :class:`~regions.RegionMeta` and
:class:`~regions.RegionVisual` for meta and visual attributes respectively.
They are subclasses of the python dictionary (`~dict`).

The meta attribute provides additional data about regions such as labels, tags,
comments, name, etc. which are used for non-display tasks.
It also stores the spectral dimensions of the region.
These classes, for now, just check whether the key is valid or not.

The valid keys for :class:`~regions.RegionMeta` class are:

1. ``label``:

    - CRTF, DS9 (text label for a region)

    - Ex: meta['label'] = 'this is a circle'

2. ``tag``:

    - DS9 (All regions may have zero or more tags associated with it,
          which may be used for grouping and searching.)

    - Ex: meta['tags'] = ['{Group 1}', '{Group 2}']}

3. ``include``:

    - CRTF, DS9 (Region inclusion)

    - Possible Values: True, False

    - Ex: meta['include'] = True

4. ``frame``:

    - CRTF (Frequency/Velocity Axis)

    - Possible values: 'REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO', 'GALACTO', 'LGROUP', 'CMB'

    - Default: image value

    - Ex: meta['frame'] = 'TOPO'

5. ``range``:

    - CRTF (Frequency/Velocity Range)

    - Possible units: GHz, MHz, kHz, km/s, Hz, channel, chan (=channel)

    - Default: image range

    - Format: [min, max]

    - Ex: meta['range'] = [-320 * u.m/u.sec, -330 * u.m/u.s]

6. ``veltype``:

    - CRTF (Velocity Calculation)

    - Possible values: 'RADIO', 'OPTICAL', 'Z', 'BETA', 'GAMMA'

    - Default: image value

    - Ex: meta['veltype'] = 'RADIO'

7. ``restfreq``:

    - CRTF (Rest Frequency)

    - Possible values: `~astropy.units.Quantity` object

    - Default: image value

    - Ex: meta['restfreq'] = Quantity("1.42GHz")

8. ``corr``:

    - CRTF (Correlational Axis)

    - Possible values: 'I', 'Q', 'U', 'V', 'RR', 'RL', 'LR', 'LL', 'XX', 'XY',
      'YX', 'YY', 'RX', 'RY', 'LX', 'LY', 'XR', 'XL', 'YR', 'YL', 'PP', 'PQ',
      'QP', 'QQ', 'RCircular', 'LCircular', 'Linear', 'Ptotal', 'Plinear',
      'PFtotal', 'PFlinear', 'Pangle'

    - Default: all planes present in image

    - Ex: meta['corr'] = ['X', 'Y']

9. ``comment``:

    - DS9, CRTF (Comment on the region)

    - Ex: meta['comment'] = 'Any comment for the region'

11. ``line``:

    - DS9 (The line region may be rendered with arrows, one at each end.
                   To indicate arrows, use the line property. A '1' indicates an arrow, '0'
                   indicates no arrow.)

    - Ex: meta['line'] = [1, 1]

12. ``name``

13. ``select``

14. ``highlite``

15. ``fixed``

16. ``edit``:

    - DS9 (The Edit property specifies if the user is allowed to edit
                  the region via the GUI.)

    - Ex: meta['edit'] = 1

17. ``move``:

    - DS9 (The Move property specifies if the user is allowed to move
          the region via the GUI. )

    - Ex: meta['move'] = 1

18. ``rotate``:

    - DS9 (The Rotate property specifies if the user is allowed to
                     rotate the region via the GUI. )

    - Ex: meta['rotate'] = 1

19. ``delete``:

    - DS9 (The Delete property specifies if the user is allowed to
                     delete the region via the GUI. )

    - Ex: meta['delete'] = 1

20. ``source``

21. ``background``



The visual attributes are meta data meant to be used to visualize regions, especially
used by plotting libraries such as `Matplotlib`_ .

The valid keys for `~regions.RegionVisual` class are:

1. ``color``: CRTF, DS9 (Region, symbol and text color)
    - Possible values: any color recognized by `Matplotlib`_, including hex values
    - Default: color=green
    - Ex: visual['color'] = 'blue'

2. ``dash``: Render region using dashed lines using current dashlist value.

3. ``font``: Name of the font.

4. ``dashlist``: Sets dashed line parameters. This does not render the region in dashed lines.

5. ``symsize``: Size of the symbol

6. ``symthick``: Thickness of the symbol

7. ``fontsize``: Size of the font.

8. ``fontstyle``: Style of the font.

9. ``usetex``: Boolean value whether the label uses tex.

10. ``labelpos``: position of the label

11. ``labeloff``: label offset

12. ``linewidth``: width of the line

13. ``linestyle``: style of the line

14. ``fill``: Boolean value whether the regions is filled

15. ``line``: The line region may be rendered with arrows, one at each end.
              To indicate arrows, use the line property. A '1' indicates an arrow, '0' indicates no arrow.

16. ``symbol``/``point``: CRTF, DS9 (Symbol for which a point region is described)

    - Ex: meta['symbol'] = 'point marker'

.. _sh-lists:

Lists
-----

A `~regions.Region` object can only represent one region, not
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
