.. _gs-contain:

Checking for points inside regions
==================================

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

    >>> from regions import PixCoord
    >>> PixCoord(55, 40) in pix_reg
    True
    >>> PixCoord(55, 200) in pix_reg
    False

The ``in`` operator only works for scalar coordinates, because Python requires
the return value to be a scalar bool, and only works for pixel regions. If you
try to use ``in`` for non-scalar coordinates, you'll get a ``ValueError``:

.. code-block:: python

    >>> pixcoord = PixCoord([50, 50], [10, 60])
    >>> pixcoord in pix_reg
    ValueError: coord must be scalar. coord=<PixCoord [(50.0, 10.0), (50.0, 60.0)]>

If you have arrays of coordinates, use the `regions.SkyRegion.contains` or
`regions.PixelRegion.contains` methods:

.. code-block:: python

    >>> pixcoords = PixCoord.from_sky(skycoords, wcs)
    >>> pix_reg.contains(pixcoords)
    array([ True, False], dtype=bool)

Note that `regions.SkyRegion.contains`
requires a WCS to be passed:

.. code-block:: python

    >>> skycoord = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sky_reg.contains(skycoord, wcs)
    array([ True, False], dtype=bool)
