0.11 (unreleased)
=================

General
-------

- The minimum required Python is now 3.11. [#598]

- The minimum required NumPy is now 1.25. [#598]

- The minimum required Astropy is now 6.0. [#598]

- The minimum required Matplotlib is now 3.8. [#598]

- The minimum required Shapely is now 2.0. [#601]

New Features
------------

Bug Fixes
---------

API Changes
-----------


0.10 (2024-09-27)
=================

General
-------

- Added support for NumPy 2.1.

API Changes
-----------

- The ``regions.test`` function has been removed. Instead use the
  ``pytest --pyarg regions`` command. [#1725]


0.9 (2024-04-12)
================

General
-------

- Added support for NumPy 2.0.

- The minimum required Python is now 3.10. [#547]

- The minimum required NumPy is now 1.23. [#547]

- The minimum required Astropy is now 5.1. [#547]

- The minimum required pytest-astropy is now 0.11. [#547]

- The minimum required sphinx-astropy is now 1.9. [#547]

Bug Fixes
---------

- Fixed a bug where the text string of ``TextPixelRegion`` and
  ``TextSkyRegion`` was dropped when serializing or writing
  to a DS9 file. [#548]


0.8 (2023-11-16)
================

General
-------

- The minimum required Python is now 3.9. [#517]

- The minimum required NumPy is now 1.22. [#517]

- The minimum required Matplotlib is now 3.5. [#517]

New Features
------------

- The ``Regions`` class can now be initialized without any arguments.
  [#527]

- The ``Regions`` ``extend`` method now can accept another ``Regions``
  object as input. [#527]

Bug Fixes
---------

API Changes
-----------

- The ``Regions`` class and its ``append`` and ``extend`` methods now
  raise a ``TypeError`` for invalid inputs. [#527]


0.7 (2022-10-27)
================

New Features
------------

- Region meta and visual metadata can be now input as a ``dict``. [#462]

Bug Fixes
---------

- Fixed an issue with the CRTF for certain malformed input. [#448]

- Fixed a bug where DS9 serialization of an empty list raised an error.
  [#449]

- Added validation to ``RegionMeta`` and ``RegionVisual`` ``update`` and
  ``setdefault`` methods. [#463]

- Fixed an issue reading CRTF file with labelcolor defined. [#473]

- Fixed a bug parsing CRTF files with ":" and "." coordinate separators.
  [#407]

API Changes
-----------

- Removed the deprecated ``BoundingBox`` class. Use
  ``RegionBoundingBox`` instead. [#456]


0.6 (2022-03-21)
================

New Features
------------

- Added the DS9 'boxcircle' point symbol. [#387]

- Added the ability to add and subtract ``PixCoord`` objects. [#396]

- Added an ``origin`` keyword to ``PolygonPixelRegion`` to allow
  specifying the vertices relative to an origin pixel. [#397]

- Added a ``RegularPolygonPixelRegion`` class. [#398]

- Added the ability to compare ``Region`` objects for equality. [#421]

- Added a ``copy`` method to ``RegionMeta`` and ``RegionVisual``. [#424]

- Added the ability to handle DS9 elliptical and rectangular annuli and
  multi-annuli regions. [#436]

- Added the ability to handle DS9 composite regions. [#436]

- Added the ability to serialize sky regions in the DS9 file format using
  the coordinate frame of the ``SkyCoord`` object. [#436]

- Added the ability to serialize multiple regions in the DS9 file format
  that have different coordinate frames. [#436]

- Added the ability to parse FITS region ``rectangle`` and
  ``rotrectangle`` shapes. [#444]

Bug Fixes
---------

- Fixed the DS9 default point symbol to use 'boxcircle'. [#387]

- Point symbol markers are no longer filled for consistency with DS9.
  [#387]

- Fixed an issue where plotting elliptical regions were incorrectly
  filled by default. [#389]

- Fixed an issue where compound region colors were being set correctly.
  [#389]

- Fixed an issue where sky/pixel conversions did not preserve ``meta``
  and ``visual`` data for some regions. [#420, #424]

- Fixed an issue with elliptical and rectangular annulus regions where
  the outer width/height could be smaller than the inner width/height.
  [#425]

- Fixed an issue converting elliptical and rectangular annulus regions
  between pixel and sky regions. [#425]

- Fixed the string representations of ``TextPixelRegion`` and
  ``TextSkyRegion`` to include quotes around the text parameter value.
  [#429]

- Fixed many issues with the DS9 parser and serializer in not
  consistently handling or preserving the region coordinate frame
  or region parameter units. [#436]

- Fixed handling of FITS shapes that are preceded by an exclamation
  mark. [#444]

- Fixed a bug where written FITS region files could not be read back in.
  [#444]

API Changes
-----------

- Removed the following deprecated I/O classes and functions:
  ``crtf_objects_to_string``, ``ds9_objects_to_string``,
  ``fits_region_object_to_table``, ``read_crtf``, ``read_ds9``,
  ``read_fits``, ``read_fits_region``, ``write_crtf``, ``write_ds9``,
  ``write_fits``, ``write_fits_region`` ``CRTFParser``, ``DS9Parser``,
  ``FITSRegionParser``, ``ShapeList``, and ``Shape``. The ``Regions``
  and ``Region`` objects now support this functionality via a unified
  I/O interface. [#386]

- Removed the deprecated ``BoundingBox`` ``slices`` attribute. [#386]

- The default matplotlib keywords that are used when plotting now depend
  on the value of ``self.visual['default_style']``. This keyword is
  currently set (to a value of 'ds9') only when reading DS9 region
  files. If set to 'ds9', DS9 plotting defaults are used. If not set or
  set to 'mpl' or None, then the matplotlib defaults will be used, with
  the exception that fill is turned off for Patch and Line2D artists.

- Renamed the ``BoundingBox`` class to ``RegionBoundingBox``. The old
  name is deprecated. [#427]

- A ``ValueError`` is raised if the radius, width, or height region
  parameters are not strictly positive (> 0). [#430]

- Added a ``precision`` keyword to the DS9 serializer and writer to
  specify the number of decimal places in output numbers. [#436]

- The ``errors`` keyword was removed from the DS9 parser and reader and
  the ``coordsys``, ``radunit``, and ``fmt`` keywords were removed from
  the DS9 serializer and writer.  The new ``precision`` keyword can be
  used when serializing and writing DS9 regions. [#436]

- The ``PixelRegion.plot()`` method now returns a
  ``matplotlib.artist.Artist`` object, which can be used in plot legends.
  [#441]

- FITS region files are now always parsed and serialized as
  ``PixelRegion`` objects. They can be converted to ``SkyRegion``
  objects using a WCS object. [#444]


0.5 (2021-07-20)
================

General
-------

- The infrastructure of the package has been updated in line with the
  APE 17 guidelines. The main changes are that the ``python setup.py
  test`` and ``python setup.py build_docs`` commands will no longer
  work. The easiest way to replicate these commands is to install the
  tox package and run ``tox -e test`` and ``tox -e build_docs``. It is
  also possible to run pytest and sphinx directly. Other significant
  changes include switching to setuptools_scm to manage the version
  number, and adding a ``pyproject.toml`` to opt in to isolated builds
  as described in PEP 517/518. [#315]

- Bump the minimum required version of Astropy to 3.2.

New Features
------------

- Added a ``as_mpl_selector`` method to the rectangular and ellipse
  pixel-based regions. This method returns an interactive Matplotlib
  selector widget. [#317]

- Added a ``get_overlap_slices`` method to ``BoundingBox``. [#348]

- Added a ``center`` attribute to ``BoundingBox``. [#348]

- Added ``get_overlap_slices`` method to ``RegionMask``. [#350]

- Added ``get_values`` method to ``RegionMask``. [#351, #353]

- Added a ``Regions`` class with a unified I/O interface for reading,
  writing, parsing, and serializing regions. [#378]

- Added ``serialize`` and ``write`` methods to all ``Region``
  subclasses. [#378]

Bug Fixes
---------

- Fixed an issue where ``RegionMask.multiply`` ``fill_value`` was not
  applied to pixels outside of the mask, but within the region bounding
  box. [#346]

- Fixed an issue where ``RegionMask.cutout`` would raise an error if
  ``fill_value`` was non-finite and the input array was integer type.
  [#346]

- A ``ValueError`` is now raised when calling ``BoundingBox.slices``
  when ``ixmin`` or ``iymin`` is negative. [#347]

- Fixed an issue in the DS9 parser where uppercase coordinate frames
  would fail. [#237]

- Fixed an issue where the CRTF file parser would fail if the CRTF
  version number was included on the first line. [#240]

- Fixed an issue where the CRTF file parser would fail if there was a
  space after the region name. [#271]

- Fixed an issue where the CRTF file parser was too restrictive about
  requiring the last and first polynomial coordinates to be the same.
  [#359, #362]

- Fixed a bug where an ``EllipsePixelRegion`` with zero height and/or
  width would raise a ``ValueError`` when creating a ``RegionMask``.
  [#363]

- Fixed parsing CRTF regions files that do not have a comma after the
  region. [#364]

- Fixed parsing CRTF regions files that contain a ``symthick`` value.
  [#365]

- Fixed an issue where ``PointPixelRegion`` objects would not plot.
  [#366]

- Fixed an issue where DS9 annulus regions with more than one annulus
  would not be parsed correctly. Such regions are skipped for now. [#371]

- Fixed an issue where ``Angle`` values for ``SkyRegion`` shape
  parameters could be incorrectly serialized. [#380]

- Fixed an issue where a semicolon in the DS9 text field would raise an
  error. [#381,#383]

- Fixed an issue where DS9 regions without metadata would not be parsed
  correctly. [#382]

- Fixed an issue parsing spaces in DS9 region metadata. [#384]

API Changes
-----------

- Deprecated the ``BoundingBox`` ``slices`` attribute. [#348]

- The ``RegionMeta`` and ``RegionVisual`` classes have been moved to the
  ``regions.core.metadata`` module. [#371]

- Deprecated the ``read_fits_region`` and ``write_fits_region``
  functions. Instead, use the ``read_fits`` and ``write_fits``
  functions. Note that the ``write_fits`` function is called as
  ``write_fits(regions, filename)`` for consistency with the other
  functions that write files. [#376]

- The following helper functions were removed from the public API:
  ``to_shape_list``, ``to_crtf_meta``, ``to_ds9_meta``,
  ``CRTFRegionParser``, ``DS9RegionParser``, ``CoordinateParser``,
  and ``FITSRegionRowParser``. [#375]

- Deprecated the following I/O classes and functions:
  ``crtf_objects_to_string``, ``ds9_objects_to_string``,
  ``fits_region_object_to_table``, ``read_crtf``, ``read_ds9``,
  ``read_fits``, ``write_crtf``, ``write_ds9``, ``write_fits``,
  ``CRTFParser``, ``DS9Parser``, ``FITSRegionParser``, ``ShapeList``,
  and ``Shape``. The ``Regions`` and ``Region`` objects now support this
  functionality via a unified I/O interface. [#378]

- Existing ``ds9`` and ``crtf`` region files will not be overwritten
  by default with the ``write`` functions. Set ``overwrite=True`` to
  overwrite existing files. [#378]


0.4 (2019-06-17)
================

New Features
------------

- Add region copy methods [#269]
- Add pixel region rotate method [#265]
- Added ``union`` and ``intersection`` methods to the ``BoundingBox``
  class. [#277]
- Add support for BOX in FITS regions [#255]
- Add PixCoord.xy [#247]

Bug Fixes
---------

- Fixed a corner-case issue where ``RegionMask.multiply()`` would not set
  non-finite data values outside of the mask but within the bounding box
  to zero. [#278]
- Fix 'text' renamed to 'label' [#234]

Other
-----

- Remove astropy-healpix dependency [#258]
- Use standalone six to avoid deprecation warnings [#235]
- Change CRTF writer to match CASA implementation [#226]
- Simplify annulus regions [#279]

See also: `regions v0.4 merged pull requests list on Github <https://github.com/astropy/regions/pulls?q=is%3Apr+milestone%3A0.4+>`__.


0.3 (2018-09-09)
================

New Features
------------

- Changed ``as_patch`` to ``as_artist`` to accommodate non-patch artists [#218]

- Implemented ``to_pixel`` for ``regions.CompoundSkyRegions``,
  ``to_mask`` for ``regions.CompoundPixelRegion`` and ``to_pixel`` for
  ``regions.CircleSkyRegion``. [#137]

- Handling dimension and broadcast of ``x`` and ``y`` in
  ``regions.PixCoord``. [#172]

- Deserialization of ``CRTF`` file format is possible. [#173]

- Added ``regions.TextPixelRegion`` and ``regions.TextSkyRegion``. [#177]

- Added ``Shape`` layer in the serialization of ``DS9`` format. Also,
  implemented ``RegionMeta`` and ``RegionVisual`` to validate
  the meta parameters. [#179]

- Serialization of ``regions.Region`` object to ``CRTF`` format
  is possible. [#186]

- Fix mask bug for regions with negative indices. [#190]

- Improved the ``plot`` methods for several regions. Added ``as_patch`` for
  annulus regions. Now, uses the parameters in the ``visual`` attributes of
  regions in the matplotlib plotting. Also, added ``mpl_properties_default``
  method in ``regions.PixelRegion`` to set the visual parameters to that of
  ``DS9`` by default. [#194]

- Now, ``to_mask`` in ``regions.CompoundPixelRegion`` handles negative
  bounding box. [#195]

- Added ``regions.RectangleAnnulusPixelRegion``,
  ``regions.RectangleAnnulusSkyRegion``, ``regions.EllipseAnnulusPixelRegion``
  and ``regions.RectangleAnnulusSkyRegion``. Also, implemented custom descriptor
  classes for attribute validation. [#196]

- Implemented FITS Region Binary Table reader and writer. [#198]

- Renamed ``Mask`` class to ``RegionMask`` and added ``origin`` arg to
  ``as_patch`` and ``plot`` methods in ``regions.Region`` class. [#203]

- Support for explicit formatting directives in ``DS9``. [#204]

See also: `regions v0.3 merged pull requests list on Github <https://github.com/astropy/regions/pulls?q=is%3Apr+milestone%3A0.3+>`__.


0.2 (2017-02-16)
================

Changelog wasn't filled.

See also: `regions v0.2 merged pull requests list on Github <https://github.com/astropy/regions/pulls?q=is%3Apr+milestone%3A0.2+>`__.


0.1 (2016-07-26)
================

Changelog wasn't filled.

See also: `regions v0.1 merged pull requests list on Github <https://github.com/astropy/regions/pulls?q=is%3Apr+milestone%3A0.1+>`__.
