
DS9 Region File Format Limitations
==================================

Region Shapes
-------------

The following DS9 regions shapes do not have corresponding
`~regions.Region` classes and therefore are not supported:

  * panda, epanda, and bpanda
  * vector
  * ruler
  * compass
  * projection

Coordinate Frames
-----------------

* `~regions.Region` objects represent abstract shapes that are not tied
  to any particular image or WCS transformation. Therefore, the following
  DS9 coordinate frames are not supported:

    * physical
    * detector
    * linear
    * amplifier
    * tile
    * wcs and wcs0
    * wcs[a-z]

  If you have a DS9 region file that uses one of these coordinate
  frames, you can load it into DS9 with its corresponding image/WCS and
  then save the region with a coordinate frame supported by this package
  (e.g., image, ICRS, FK5, etc.).

* A coordinate frame must be specified in the DS9 region file. This
  package does not assume a default coordinate frame.

* DS9 regions with mixed coordinate systems within a region specifier
  are not supported. For example, a region file with an ``image``
  coordinate frame must not have angular sizes, e.g., ``image; circle(650,
  932, 3')`` is not supported.


Other Limitations
-----------------

* Reading and then writing a DS9 region file will not produce an
  identical file to the original, but the encoded regions are identical.
  Therefore, it will produce identical `~regions.Region` objects
  when read back in again. In other words, read/write/read (or
  parse/serialize/parse) will exactly roundtrip `~regions.Region`
  objects.

* DS9 composite regions are parsed into separate, independent
  `~regions.Region` objects.


Plotting Differences
--------------------

* The point symbol size must be a integer (matplotlib limitation).

* Point symbols may not have dashed lines (matplotlib limitation).

* Point symbols can be plotted with filled colors, which DS9 does not
  support. Because of this, DS9 serializations will ignore the
  ``fill=1`` metadata.

* Annulus regions can be plotted with filled colors, which DS9 does not
  support. Because of this, DS9 serializations will ignore the
  ``fill=1`` metadata for annulus regions (annulus, ellipse, box).

* Text labels are plotted only for `~regions.TextPixelRegion` objects.
  For other region objects, the text labels are stored in the region
  metadata.
