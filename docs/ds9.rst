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
    Type: reg
    Coord sys: galactic
    Region type: circle
    Meta: {'color': 'green', 'include': True}
    Composite: False
    Include: True
    >>> regions = parser.shapes.to_regions()
    >>> print(regions[0])
    Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (42., 43.)>
    radius: 3.0 deg

Serialisation is done using the `~regions.ds9_objects_to_string` function

    >>> from regions import ds9_objects_to_string
    >>> ds9_objects_to_string(regions, coordsys='galactic')
    '# Region file format: DS9 astropy/regions\ngalactic\ncircle(42.000000,43.000000,3.000000) # color=green\n'

There's also `~regions.write_ds9` and `~regions.read_ds9` which write to and
read from a file in addition to doing the region serialisation and parsing.

.. code-block:: python

    >>> from regions import read_ds9, write_ds9
    >>> filename = 'ds9.reg'
    >>> write_ds9(regions, filename)
    >>> regions = read_ds9(filename)
    >>> regions
    [<CircleSkyRegion(center=<SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        (245.347655, 24.429081)>, radius=3.0 deg)>]

The ``visual`` metadata includes items used for display, e.g.:

.. code-block:: python

    >>> print(regions[0].visual)
    {'color': 'green'}

Some of these keyword may eventually be used by the plotting utilities and
standardized, but they are not as of v0.3.
