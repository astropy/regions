.. _gs-ds9:

Reading/writing to DS9 region files
===================================

The regions package provides the functionality to serialise and de-serialise
Python lists of `~regions.Region` objects to DS9 region strings.
De-serialisation is done using  the `~regions.DS9Parser`. It converts the DS9
string to `~regions.ShapeList` object, which is a list of `~regions.Shape` each
representing one DS9 region. The `~regions.Shape` objects can be converted to
`~regions.Region` objects.

.. code-block:: python

    >>> from regions import DS9Parser
    >>> reg_string = 'galactic\ncircle(42,43,3) # color=green'
    >>> parser = DS9Parser(reg_string)
    >>> print(parser.shapes[0])
    Shape
    Format Type : DS9
    Coord sys : galactic
    Region type : circle
    Coord: [<Angle 42.0 deg>, <Angle 43.0 deg>, <Quantity 3.0 deg>]
    Meta: {'color': 'green', 'include': ''}
    Composite: False
    Include: True
    >>> regions = parser.shapes.to_regions()
    >>> print(regions[0])
    Region: CircleSkyRegion
    center: <Galactic Coordinate: (l, b) in deg
        ( 42.,  43.)>
    radius: 3.0 deg

Serialisation is done using the `~regions.ds9_objects_to_string` function

    >>> from regions import ds9_objects_to_string
    >>> ds9_objects_to_string(regions, coordsys='galactic')
    '# Region file format: DS9 astropy/regions\ngalactic\ncircle(42.0000,43.0000,3.0000)\n'

.. warning::

    The API for serialisation to DS9 string is currently not using the
    `~regions.Shape` layer. It should be adapted.  Also, all regions currently
    have ``meta`` and ``visual`` arguments for ``__init__`` and stored as
    region data members. These need to be documented and tests added, or
    removed.

There's also `~regions.write_ds9` and `~regions.read_ds9` with write to and
read from a file in addition to doing the region serialisation and parsing.

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

The ``visual`` metadata includes items used for display, e.g.::

    >>> r.visual
    {'color': 'green',
     'font': '"helvetica 10 normal roman" ',
     'point': 'x',
     'width': '1'}

Some of these keyword may eventually be used by the plotting utilities and
standardized, but they are not as of v0.3.
