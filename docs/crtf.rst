.. _gs-crtf:

Reading/writing to CRTF region files
====================================

CASA Region Text Format (CRTF) : The CASA region file format provides a
flexible, easily edited set of region definitions which are accepted across
CASA tasks.

The regions package provides the functionality to serialise and de-serialise
Python lists of `~regions.Region` objects to CRTF region strings.
De-serialisation is done using  the `~regions.CRTFParser`. It converts the CRTF
string to `~regions.ShapeList` object, which is a list of `~regions.Shape` each
representing one CRTF region. The `~regions.Shape` objects can be converted to
`~regions.Region` objects.

To learn more about the format specifications, please go to `the official page
<https://casaguides.nrao.edu/index.php/CASA_Region_Format#Region_definitions>`_.

.. code-block:: python

    >>> from regions import CRTFParser
    >>> reg_string = 'circle[[42deg, 43deg], 3deg], coord=J2000, color=green '
    >>> parser = CRTFParser(reg_string)
    >>> print(parser.shapes[0])
    Shape
    Type : reg
    Coord sys : fk5
    Region type : circle
    Meta: {'color': 'green', 'include': True, 'type': 'reg'}
    Composite: False
    Include: True
    >>> regions = parser.shapes.to_regions()
    >>> print(regions[0])
    Region: CircleSkyRegion
    center: <SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        (42., 43.)>
    radius: 3.0 deg
    >>> print(regions[0].meta)
    {'include': True, 'type': 'reg'}
    >>> print(regions[0].visual)
    {'color': 'green'}

Serialisation is done using the `~regions.crtf_objects_to_string` function

.. code-block:: python

    >>> from regions import crtf_objects_to_string
    >>> crtf_objects_to_string(regions, coordsys='galactic')
    '#CRTFv0\nglobal coord=galactic\ncircle[[144.559169deg, -14.923593deg], 3.000000deg], color=green\n'

There's also `~regions.write_crtf` and `~regions.read_crtf` which write to and
read from a file in addition to doing the region serialisation and parsing.

.. code-block:: python

    >>> from regions import read_crtf, write_crtf
    >>> filename = 'region.crtf'
    >>> write_crtf(regions, filename)
    >>> regions = read_crtf(filename)
    >>> regions
    [<CircleSkyRegion(center=<SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
        (42., 43.)>, radius=3.0 deg)>]
