.. _stcs_io:

STC-S Region File Format
========================

The STC-S (Space-Time Coordinate String) format is an IVOA (International
Virtual Observatory Alliance) standard for describing spatial and temporal
regions and coordinates in a string representation. This format provides a
compact, human-readable way to express astronomical regions.

Introduction
------------

STC-S is defined in the `IVOA STC-S Note`_ and provides a standardized way to
represent spatial regions. The format is particularly useful for:

* Interoperability between astronomical software
* Compact representation of regions in databases
* Human-readable region definitions
* Exchange of regions between observatories and data centers

The regions package supports reading from and writing to STC-S files and
parsing/serializing STC-S strings.

Features
--------

- **Reading**: Parse STC-S strings and files into `regions.Region` objects
- **Writing**: Serialize `regions.Region` objects to STC-S format
- **Round-trip**: Full round-trip conversion preserving region properties
- **Multiple Shapes**: Support for circles, ellipses, boxes/rectangles, polygons, and points
- **Coordinate Systems**: Support for ICRS, FK5, FK4, Galactic, Ecliptic, and Image coordinates
- **Reference Positions**: Support for various reference positions (Barycenter, Geocenter, Topocenter, etc.)

Usage Examples
--------------

Reading STC-S Files
^^^^^^^^^^^^^^^^^^^

To read an STC-S region file, you **must** specify the format explicitly,
as STC-S files cannot be reliably auto-detected:

.. doctest-skip::

    >>> from regions import Regions
    >>> regions = Regions.read('my_regions.stcs', format='stcs')

.. note::
    Unlike DS9 and CRTF formats, STC-S files do not have unique file
    signatures and cannot be auto-detected based on content. You must
    always specify ``format='stcs'`` when reading STC-S files.

Parsing STC-S Strings
^^^^^^^^^^^^^^^^^^^^^

You can also parse STC-S strings directly:

.. doctest-skip::

    >>> stcs_string = """
    ... # Example STC-S regions
    ... # Examples created for demonstration based on the IVOA STC-S standard
    ... Circle ICRS BARYCENTER 180.0 10.0 0.5
    ... Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
    ... Position FK5 GEOCENTER 85.0 -15.0
    ... """
    >>> regions = Regions.parse(stcs_string, format='stcs')
    >>> print(regions)
    [<CircleSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in deg (180., 10.)>, radius=0.5 deg)>,
     <EllipseSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in deg (150., -20.)>, width=2.0 deg, height=1.0 deg, angle=45.0 deg)>,
     <PointSkyRegion(center=<SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg (85., -15.)>)>]

Writing STC-S Files
^^^^^^^^^^^^^^^^^^^

To write regions to an STC-S file:

.. doctest-skip::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from regions import CircleSkyRegion, Regions
    >>> center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
    >>> circle = CircleSkyRegion(center=center, radius=0.5 * u.degree)
    >>> regions = Regions([circle])
    >>> regions.write('output.stcs', format='stcs', overwrite=True)

Serializing to STC-S Strings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To serialize regions to an STC-S string:

.. doctest-skip::

    >>> stcs_string = regions.serialize(format='stcs')
    >>> print(stcs_string)
    Circle ICRS BARYCENTER 180 10 0.5

Working with Pixel Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

STC-S also supports pixel coordinates using the IMAGE frame:

.. doctest-skip::

    >>> from regions.core import PixCoord
    >>> from regions import CirclePixelRegion
    >>> center = PixCoord(100.5, 200.3)
    >>> pixel_circle = CirclePixelRegion(center=center, radius=15.0)
    >>> stcs_string = pixel_circle.serialize(format='stcs')
    >>> print(stcs_string)
    Circle IMAGE UNKNOWN 100.5 200.3 15

Complete Examples
^^^^^^^^^^^^^^^^^

Reading STC-S Files and Strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doctest-skip::

    >>> from regions import Regions
    >>>
    >>> # Parse STC-S string directly
    >>> stcs_string = \"\"\"
    ... # Example STC-S regions
    ... # Examples created for demonstration based on the IVOA STC-S standard
    ... Circle ICRS BARYCENTER 180.0 10.0 0.5
    ... Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
    ... Position FK5 GEOCENTER 85.0 -15.0
    ... \"\"\"
    >>>
    >>> regions = Regions.parse(stcs_string, format='stcs')
    >>> print(f"Parsed {len(regions)} regions:")
    Parsed 3 regions:
    >>> for i, region in enumerate(regions):
    ...     print(f"  {i+1}. {region}")
      1. <CircleSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in deg (180., 10.)>, radius=0.5 deg)>
      2. <EllipseSkyRegion(center=<SkyCoord (ICRS): (ra, dec) in deg (150., -20.)>, width=2.0 deg, height=1.0 deg, angle=45.0 deg)>
      3. <PointSkyRegion(center=<SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg (85., -15.)>)>

Writing Multiple Region Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doctest-skip::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from regions.core import PixCoord, Regions
    >>> from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
    ...                            EllipseSkyRegion, PolygonSkyRegion,
    ...                            PointSkyRegion)
    >>>
    >>> # Create some example regions
    >>> regions = []
    >>>
    >>> # Sky circle
    >>> center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
    >>> circle = CircleSkyRegion(center=center, radius=0.5 * u.degree)
    >>> regions.append(circle)
    >>>
    >>> # Sky ellipse
    >>> center = SkyCoord(150.0, -20.0, unit='degree', frame='icrs')
    >>> ellipse = EllipseSkyRegion(center=center, width=2.0 * u.degree,
    ...                          height=1.0 * u.degree, angle=45.0 * u.degree)
    >>> regions.append(ellipse)
    >>>
    >>> # Sky polygon
    >>> vertices = SkyCoord([45.0, 50.0, 50.0, 45.0],
    ...                    [45.0, 45.0, 50.0, 50.0],
    ...                    unit='degree', frame='icrs')
    >>> polygon = PolygonSkyRegion(vertices=vertices)
    >>> regions.append(polygon)
    >>>
    >>> # Point region
    >>> center = SkyCoord(85.0, -15.0, unit='degree', frame='fk5')
    >>> point = PointSkyRegion(center=center)
    >>> regions.append(point)
    >>>
    >>> # Pixel circle
    >>> center = PixCoord(100.5, 200.3)
    >>> pixel_circle = CirclePixelRegion(center=center, radius=15.0)
    >>> regions.append(pixel_circle)
    >>>
    >>> regions_obj = Regions(regions)
    >>>
    >>> # Serialize to STC-S string
    >>> stcs_string = regions_obj.serialize(format='stcs')
    >>> print("Serialized STC-S:")
    Serialized STC-S:
    >>> print(stcs_string)
    Circle ICRS BARYCENTER 180 10 0.5
    Ellipse ICRS BARYCENTER 150 -20 1 0.5 45
    Polygon ICRS BARYCENTER 45 45 50 45 50 50 45 50
    Position FK5 BARYCENTER 85 -15
    Circle IMAGE UNKNOWN 100.5 200.3 15
    >>>
    >>> # Write to file (uncomment to actually write)
    >>> # regions_obj.write('output.stcs', format='stcs', overwrite=True)

Round-trip Conversion
~~~~~~~~~~~~~~~~~~~~~

Verify that regions can be converted to STC-S and back without loss:

.. doctest-skip::

    >>> # Original STC-S
    >>> original_stcs = \"\"\"Circle ICRS BARYCENTER 180.0 10.0 0.5
    ... Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
    ... Position FK5 GEOCENTER 85.0 -15.0\"\"\"
    >>>
    >>> print("Original STC-S:")
    Original STC-S:
    >>> print(original_stcs)
    Circle ICRS BARYCENTER 180.0 10.0 0.5
    Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0
    Position FK5 GEOCENTER 85.0 -15.0
    >>>
    >>> # Parse
    >>> regions = Regions.parse(original_stcs, format='stcs')
    >>> print(f"\\nParsed {len(regions)} regions")

    Parsed 3 regions
    >>>
    >>> # Serialize back
    >>> serialized = regions.serialize(format='stcs')
    >>> print("\\nSerialized back to STC-S:")

    Serialized back to STC-S:
    >>> print(serialized)
    Circle ICRS BARYCENTER 180 10 0.5
    Ellipse ICRS BARYCENTER 150 -20 1 0.5 45
    Position FK5 BARYCENTER 85 -15
    >>>
    >>> # Parse again to verify
    >>> regions2 = Regions.parse(serialized, format='stcs')
    >>> print(f"\\nRe-parsed {len(regions2)} regions - round-trip successful!")

    Re-parsed 3 regions - round-trip successful!

STC-S Format Specification
---------------------------

Basic Syntax
^^^^^^^^^^^^

The basic syntax for STC-S regions is:

.. code-block::

    <Shape> <Frame> <RefPos> <Coordinates> [<Parameters>]

Where:

* **Shape**: The geometric shape (Circle, Ellipse, Box, Polygon, Position)
* **Frame**: The coordinate reference frame (ICRS, FK5, FK4, GALACTIC, ECLIPTIC, IMAGE)
* **RefPos**: The reference position (BARYCENTER, GEOCENTER, TOPOCENTER, etc.)
* **Coordinates**: The shape-specific coordinate parameters
* **Parameters**: Additional shape parameters (radii, angles, etc.)

Supported Shapes
^^^^^^^^^^^^^^^^

Circle
~~~~~~

Defines a circular region:

.. code-block::

    Circle <Frame> <RefPos> <center_lon> <center_lat> <radius>

Example:

.. code-block::

    Circle ICRS BARYCENTER 180.0 10.0 0.5

Ellipse
~~~~~~~

Defines an elliptical region:

.. code-block::

    Ellipse <Frame> <RefPos> <center_lon> <center_lat> <semi_major> <semi_minor> <angle>

Example:

.. code-block::

    Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0

Box
~~~

Defines a rectangular region:

.. code-block::

    Box <Frame> <RefPos> <center_lon> <center_lat> <width> <height> <angle>

Example:

.. code-block::

    Box ICRS BARYCENTER 120.0 30.0 2.0 1.0 0.0

Polygon
~~~~~~~

Defines a polygonal region with multiple vertices:

.. code-block::

    Polygon <Frame> <RefPos> <lon1> <lat1> <lon2> <lat2> <lon3> <lat3> ...

Example:

.. code-block::

    Polygon ICRS BARYCENTER 45.0 45.0 50.0 45.0 50.0 50.0 45.0 50.0

Position
~~~~~~~~

Defines a point region:

.. code-block::

    Position <Frame> <RefPos> <lon> <lat>

Example:

.. code-block::

    Position FK5 GEOCENTER 85.0 -15.0

Coordinate Frames
^^^^^^^^^^^^^^^^^

The following coordinate reference frames are supported:

=============  ===============================================
Frame          Description
=============  ===============================================
ICRS           International Celestial Reference System
FK5            Fifth Fundamental Catalogue (J2000.0)
FK4            Fourth Fundamental Catalogue (B1950.0)
GALACTIC       Galactic coordinate system
ECLIPTIC       Ecliptic coordinate system
IMAGE          Pixel/image coordinates
=============  ===============================================

Reference Positions
^^^^^^^^^^^^^^^^^^^

The following reference positions are supported:

=================  ===============================================
Reference Position Description
=================  ===============================================
BARYCENTER         Solar system barycenter
GEOCENTER          Earth center
TOPOCENTER         Earth surface/topocentric
HELIOCENTER        Sun center
LSR                Local Standard of Rest
LSRK               Kinematic Local Standard of Rest
LSRD               Dynamic Local Standard of Rest
UNKNOWN            Unspecified reference position
=================  ===============================================

File Format
^^^^^^^^^^^

STC-S files typically use the following extensions:

* ``.stcs``
* ``.stc``
* ``.stcs.txt``
* ``.stc.txt``

Files can contain:

* Comments starting with ``#``
* Multiple regions, one per line
* Blank lines (ignored)

All STC-S files generated by astropy-regions include a standard header:

.. code-block::

    # Region file format: STC-S astropy/regions

Example STC-S file:

.. code-block::

    # Region file format: STC-S astropy/regions
    # Examples created based on the IVOA STC-S standard specification

    # Central source
    Circle ICRS BARYCENTER 180.0 10.0 0.5

    # Extended emission
    Ellipse ICRS BARYCENTER 150.0 -20.0 1.0 0.5 45.0

    # Point sources
    Position FK5 GEOCENTER 85.0 -15.0
    Position FK5 GEOCENTER 90.0 -10.0

Region Mapping
--------------

The following table shows the mapping between STC-S shapes and
astropy-regions classes:

=============  =================================  =================================
STC-S Shape    Sky Region Class                   Pixel Region Class
=============  =================================  =================================
Circle         `~regions.CircleSkyRegion`         `~regions.CirclePixelRegion`
Ellipse        `~regions.EllipseSkyRegion`        `~regions.EllipsePixelRegion`
Box            `~regions.RectangleSkyRegion`      `~regions.RectanglePixelRegion`
Polygon        `~regions.PolygonSkyRegion`        `~regions.PolygonPixelRegion`
Position       `~regions.PointSkyRegion`          `~regions.PointPixelRegion`
=============  =================================  =================================

Format Limitations
------------------

Region Shapes
^^^^^^^^^^^^^

The following STC-S features are not currently supported:

* **Time coordinates and temporal regions**:

  .. code-block::

      Time TT TOPOCENTER 2000-01-01T12:00:00 2000-01-02T12:00:00
      TimeInterval TT GEOCENTER 2000-01-01T00:00:00 2001-01-01T00:00:00

* **Spectral coordinates**:

  .. code-block::

      Spectral TOPOCENTER 1420.4 MHz
      SpectralInterval BARYCENTER 1400.0 1440.0 MHz

* **Redshift specifications**:

  .. code-block::

      RedshiftInterval BARYCENTER VELOCITY OPTICAL 200.0 2300.0 km/s
      Redshift BARYCENTER VELOCITY RADIO 0.1

* **Complex compound operations**:

  .. code-block::

      Union ICRS BARYCENTER (Circle 180 10 0.5) (Circle 190 20 0.5)
      Intersection ICRS BARYCENTER (Circle 180 10 2.0) (Box 180 10 1.0 1.0 0.0)
      Difference ICRS BARYCENTER (Circle 180 10 2.0) (Circle 180 10 0.5)

* **Unit specifications and mixed units**:

  .. code-block::

      Circle ICRS BARYCENTER unit deg arcsec 180.0 10.0 30.0

* **Error bounds and uncertainties**:

  .. code-block::

      Circle ICRS BARYCENTER 180.0 10.0 0.5 Error 0.1 0.1 0.05

* **Resolution and pixel size specifications**:

  .. code-block::

      Circle ICRS BARYCENTER 180.0 10.0 0.5 Resolution 0.1 PixSize 0.05

Coordinate Systems
^^^^^^^^^^^^^^^^^^

* Only spatial coordinates are supported; temporal coordinates are ignored
* Complex coordinate transformations are not implemented
* Some specialized coordinate systems may not be fully supported

Auto-detection Limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^

* **STC-S files cannot be auto-detected** based on content, as they lack unique
  file signatures and use keywords that also appear in DS9, CRTF, and other
  region formats. You must always specify ``format='stcs'`` explicitly when
  reading STC-S files.

* Auto-detection only works based on file extensions (``.stcs``, ``.stc``,
  ``.stcs.txt``, ``.stc.txt``).

Other Limitations
^^^^^^^^^^^^^^^^^

* Reading and writing an STC-S file will not produce an identical file to the
  original, but the encoded regions are identical. The regions will produce
  identical `~regions.Region` objects when read back in again.

* Comments and formatting may not be preserved exactly during round-trip
  operations.

* Error handling for malformed STC-S strings could be more detailed.

Examples
--------

Complete Example
^^^^^^^^^^^^^^^^

Here's a complete example showing how to work with STC-S files:

.. doctest-skip::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from regions import (CircleSkyRegion, EllipseSkyRegion,
    ...                      PolygonSkyRegion, PointSkyRegion, Regions)

    >>> # Create some regions
    >>> regions = []

    >>> # Sky circle
    >>> center = SkyCoord(180.0, 10.0, unit='degree', frame='icrs')
    >>> circle = CircleSkyRegion(center=center, radius=0.5 * u.degree)
    >>> regions.append(circle)

    >>> # Sky ellipse
    >>> center = SkyCoord(150.0, -20.0, unit='degree', frame='icrs')
    >>> ellipse = EllipseSkyRegion(center=center, width=2.0 * u.degree,
    ...                          height=1.0 * u.degree, angle=45.0 * u.degree)
    >>> regions.append(ellipse)

    >>> # Polygon
    >>> vertices = SkyCoord([45.0, 50.0, 50.0, 45.0],
    ...                    [45.0, 45.0, 50.0, 50.0],
    ...                    unit='degree', frame='icrs')
    >>> polygon = PolygonSkyRegion(vertices=vertices)
    >>> regions.append(polygon)

    >>> # Point
    >>> center = SkyCoord(85.0, -15.0, unit='degree', frame='fk5')
    >>> point = PointSkyRegion(center=center)
    >>> regions.append(point)

    >>> regions_obj = Regions(regions)

    >>> # Write to file
    >>> regions_obj.write('example.stcs', format='stcs', overwrite=True)

    >>> # Read back
    >>> read_regions = Regions.read('example.stcs', format='stcs')
    >>> print(f"Read {len(read_regions)} regions")
    Read 4 regions

Round-trip Conversion
^^^^^^^^^^^^^^^^^^^^

STC-S supports full round-trip conversion:

.. doctest-skip::

    >>> # Original STC-S
    >>> original = "Circle ICRS BARYCENTER 180.0 10.0 0.5"
    >>>
    >>> # Parse -> Serialize -> Parse
    >>> regions = Regions.parse(original, format='stcs')
    >>> serialized = regions.serialize(format='stcs')
    >>> regions2 = Regions.parse(serialized, format='stcs')
    >>>
    >>> # Verify consistency
    >>> r1, r2 = regions[0], regions2[0]
    >>> print(f"Original center: {r1.center}")
    >>> print(f"Round-trip center: {r2.center}")
    >>> print(f"Centers match: {r1.center.separation(r2.center) < 1e-10 * u.degree}")

Implementation Details
----------------------

The STC-S module consists of:

- **`core.py`**: Core parsing functions, mappings, and utilities
- **`connect.py`**: File format identification and registry
- **`read.py`**: STC-S reading and parsing functionality
- **`write.py`**: STC-S writing and serialization functionality
- **`tests/`**: Comprehensive test suite

Key Functions
^^^^^^^^^^^^^

- `validate_stcs_string()`: Validate STC-S format
- `parse_coordinate_frame()`: Extract coordinate frame and reference position
- `parse_numbers()`: Parse numeric parameters
- `_parse_stcs()`: Main parsing function
- `_serialize_stcs()`: Main serialization function

Testing
-------

The STC-S module includes a comprehensive test suite. Run tests using:

.. code-block:: bash

    # Run all STC-S tests
    pytest regions/io/stcs/tests/

    # Run specific test file
    pytest regions/io/stcs/tests/test_stcs.py

    # Run with verbose output
    pytest regions/io/stcs/tests/ -v

Test data files are located in `regions/io/stcs/tests/data/` and include:

- `stcs_basic.stcs`: Basic shape examples (created based on IVOA STC-S standard)
- `stcs_pixel.stcs`: Pixel coordinate examples (created based on IVOA STC-S standard)
- `stcs_complex.stcs`: Complex region examples (inspired by CDS STC-S Rust implementation)

Example Sources and Attribution
-------------------------------

The STC-S examples used throughout this documentation and in test files come from the following sources:

**Test Data Files:**

- **Basic Examples** (`stcs_basic.stcs`, `stcs_pixel.stcs`): Created specifically for this implementation based on the IVOA STC-S standard specification. These examples demonstrate fundamental STC-S syntax and cover the core region types and coordinate systems.

- **Complex Examples** (`stcs_complex.stcs`): Some coordinate examples inspired by the CDS STC-S Rust implementation, adapted to test various coordinate systems and reference positions. Additional examples created based on the IVOA standard.

**Documentation Examples:**

- **Tutorial Examples**: All examples in the usage sections were created specifically for demonstration purposes, following the IVOA STC-S standard syntax to illustrate proper usage patterns.

- **Format Specification Examples**: Directly based on or derived from examples in the IVOA STC-S standard documentation to ensure accuracy and compliance.

**Advanced Syntax Examples:**

- **Unsupported Features**: Examples of advanced STC-S syntax (time coordinates, spectral coordinates, compound operations) are referenced from the IVOA STC-S standard and CDS implementations to show what the standard supports beyond the current implementation scope.

References
----------

- `IVOA STC-S Standard`_
- `CDS STC-S Rust Implementation`_
- `astropy-regions Issue #21`_

.. _IVOA STC-S Note: https://www.ivoa.net/documents/Notes/STC-S/20091030/NOTE-STC-S-1.33-20091030.html
.. _IVOA STC-S Standard: https://www.ivoa.net/documents/Notes/STC-S/20091030/NOTE-STC-S-1.33-20091030.html
.. _CDS STC-S Rust Implementation: https://github.com/cds-astro/cds-stc-rust/
.. _astropy-regions Issue #21: https://github.com/astropy/regions/issues/21
