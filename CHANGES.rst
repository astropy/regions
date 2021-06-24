0.5 (unreleased)
================

General
-------

- The infrastructure of the package has been updated in line with the
  APE 17 guidelines. The main changes are that the ``python setup.py
  test`` and ``python setup.py build_docs`` commands will no longer
  work. The easiest way to replicate these commands is to install the
  tox (https://tox.readthedocs.io) package and run ``tox -e test`` and
  ``tox -e build_docs``. It is also possible to run pytest and sphinx
  directly. Other significant changes include switching to setuptools_scm
  to manage the version number, and adding a ``pyproject.toml`` to opt in
  to isolated builds as described in PEP 517/518. [#315]

- Bump the minimum required version of Astropy to 3.2.

New Features
------------

- Added a ``as_mpl_selector`` method to the rectangular and ellipse
  pixel-based regions. This method returns an interactive Matplotlib
  selector widget. [#317]

- Added a unified read-write interface for all region formats, inspired
  by and using the same infrastructure as astropy.table. [#307]

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
  requiring the last and first poly coordinates to be the same.
  [#359, #362]

- Fixed a bug where an ``EllipsePixelRegion`` with zero height and/or
  width would raise a ``ValueError`` when creating a ``RegionMask``.
  [#363]

- Fixed parsing CRTF regions files that do not have a comma after the
  region. [#364]

- Fixed parsing CRTF regions files contain a ``symthick`` value. [#365]

- Fixed an issue where ``PointPixelRegion`` objects would not plot.
  [#366]

- Fixed an issue where DS9 annulus regions with more than one annulus
  would not be parsed correctly. Such regions are skipped for now. [#371]


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

New features
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
