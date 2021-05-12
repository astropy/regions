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

API Changes
-----------

- Deprecated the ``BoundingBox`` ``slices`` attribute. [#348]


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

- Handling dimension and broadcast of `x` and `y` in ``regions.PixCoord``.
  [#172]

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
