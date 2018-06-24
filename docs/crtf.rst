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

To know more about the format specifications, please go to `the official page
<https://casaguides.nrao.edu/index.php/CASA_Region_Format#Region_definitions>`_

.. code-block:: python

    >>> from regions import CRTFParser
    >>> reg_string = 'circle[[42deg, 43deg], 3deg], coord=J2000, color=green '
    >>> parser = CRTFParser(reg_string)
    >>> print(parser.shapes[0])
    Shape
    Format Type : CRTF
    Type: reg
    Coord sys : fk5
    Region type : circle
    Coord: [<Angle 42.0 deg>, <Angle 43.0 deg>, <Quantity 3.0 deg>]
    Meta: {'color': 'green', 'coord': 'j2000', 'type': 'reg'}
    Composite: False
    Include: True

    >>> regions = parser.shapes.to_regions()
    >>> print(regions[0])
    Region: CircleSkyRegion
    center: <FK5 Coordinate (equinox=J2000.000): (ra, dec) in deg
        ( 42.,  43.)>
    radius: 3.0 deg
