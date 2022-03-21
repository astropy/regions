
FITS Region File Format Limitations
===================================

Region Shapes
-------------

* FITS regions are specified in pixel units, thus a `~regions.SkyRegion`
  cannot be serialized to the FITS region format. Regions are always
  parsed and serialized as `~regions.PixelRegion` objects when using the
  FITS region format.

* When reading a ``elliptannulus``, the first ``ROTANG`` will be used
  for both the inner and outer ellipse. The second ``ROTANG`` will be
  ignored. In other words, you cannot have an elliptical annulus where the
  inner and outer ellipse have different rotation angles.

* Shapes where the value in the SHAPE column is preceded by an exclamation
  mark (e.g., ``!circle``) will be read in as a `~regions.PixelRegion`
  object and have ``include = 0`` in their ``meta`` dictionary. When
  such objects are serialized, their shape will be prepended by ``!``.

* The FITS ``box``, ``rotbox``, ``rectangle``, and ``rotrectangle``
  shapes will all be parsed as a `~regions.RectanglePixelRegion`. In
  turn, `~regions.RectanglePixelRegion` is always serialized as a
  ``rotbox`` shape.

* The following `~regions.PixelRegion` classes do not have corresponding
  FITS shapes and therefore are not supported (a warning will be raised):

    * `~regions.RectangleAnnulusPixelRegion`
    * `~regions.LinePixelRegion`
    * `~regions.TextPixelRegion`
    * `~regions.CompoundPixelRegion`

* FITS regions are always parsed and serialized into separate regions.
  Shapes that have the ``COMPONENT`` column will have that value
  stored in the `~regions.PixelRegion` ``meta`` dictionary with the
  ``component`` key. Such regions will include the ``COMPONENT`` column
  when serialized.

* FITS parsing and serialization use only the ``include`` and
  ``component`` metadata and no visual metadata.


Coordinate Frames
-----------------

* `~regions.Region` objects represent abstract shapes that are not
  tied to any particular image or WCS transformation. Therefore, any
  WCS information in the FITS region file header will not be read.
  Regions are always parsed and serialized as `~regions.PixelRegion`
  objects when using the FITS region format. However, if desired you
  can use the WCS information in the FITS region file to convert a
  `~regions.PixelRegion` object to a `~regions.SkyRegion` object, e.g.,:

.. doctest-skip::

    >>> from astropy.io import fits
    >>> from astropy.wcs import WCS
    >>> header = fits.getheader('my_region.fits', 1)
    >>> wcs = WCS(header, keysel=('image', 'binary', 'pixel'))
    >>> sky_region = pix_region.to_sky(wcs)


Other Limitations
-----------------

* Reading and then writing a FITS region file will not produce an
  identical file to the original, but the encoded regions are identical.
  Therefore, it will produce identical `~regions.Region` objects
  when read back in again. In other words, read/write/read (or
  parse/serialize/parse) will exactly roundtrip `~regions.Region`
  objects.
