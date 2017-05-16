.. _gs-ds9:

Reading/writing to DS9 region files
===================================

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
