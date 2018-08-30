0.3 (Unreleased)
================

NEW FEATURES
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


0.2 (16-02-17)
==============


0.1 (26-07-16)
==============

    -
