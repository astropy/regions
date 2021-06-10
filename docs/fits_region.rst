.. _gs-fits:

Reading/writing to FITS region files
====================================

The regions package provides the functionality to serialise
and de-serialise Python lists of `~regions.Region`
objects to FITS region file. A `FITS Region Binary Table
<https://fits.gsfc.nasa.gov/registry/region.html>`_ defines a spatial
region of a two-dimensional image in pixels.

The file can be read in using the `~regions.read_fits_region`
function. The name of the header must be ``REGION`` for
the `~regions.read_fits_region` to parse the table. The
`~regions.read_fits_region` function returns a sky region object.

Some FITS regions, such as ``rectangle``, ``rotrectangle``,
``pie``, ``sector``, are not supported by this package.

The file can also be read into a `~astropy.table.Table` object using the
`~astropy.table.Table.read` method::

    >>> from astropy.table import Table
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('data/fits_region.fits',
    ...                                  package='regions.io.fits.tests')
    >>> table = Table.read(filename)
    >>> print(table)
         X [4]            Y [4]        SHAPE       R [4]         ROTANG  COMPONENT
          pix              pix                      pix           deg
    ---------------- ---------------- ------- ---------------- --------- ---------
       2896.5 .. 0.0    5056.5 .. 0.0  circle  381.9716 .. 0.0       0.0         1
    5282.0541 .. 0.0 4854.5699 .. 0.0  ROTBOX 1303.4597 .. 0.0 28.395178         2
       2944.5 .. 0.0    3472.5 .. 0.0 ellipse     288.0 .. 0.0  337.4048         3
        341.0 .. 0.0     345.0 .. 0.0  ROTBOX      56.0 .. 0.0      65.0         4
        341.0 .. 0.0     345.0 .. 0.0 ANNULUS      56.0 .. 0.0       0.0         5
        341.0 .. 0.0     345.0 .. 0.0   point       0.0 .. 0.0       0.0         6
        341.0 .. 0.0     345.0 .. 0.0   point       0.0 .. 0.0       0.0         7
          1.0 .. 4.0       5.0 .. 8.0 polygon       0.0 .. 0.0       0.0         8
         10.0 .. 0.0       5.5 .. 0.0     BOX      10.0 .. 0.0       0.0         9


The :class:`~regions.FITSRegionParser` class can be used to parse this
table. It converts the FITS region table to a `~regions.ShapeList`
object, which is a list of `~regions.Shape` objects, each one
representing a FITS region in pixels::

    >>> from regions import FITSRegionParser
    >>> parser = FITSRegionParser(table)
    >>> print(parser.shapes[0])
    Shape
    Type: reg
    Coord sys: physical
    Region type: circle
    Meta: {'tag': '1'}
    Composite: False
    Include: False

Finally, a list of `~regions.Region` objects can be
created from the `~regions.ShapeList` object using the
:meth:`~regions.ShapeList.to_regions` method::

    >>> regions = parser.shapes.to_regions()  # pixel regions
    >>> print(regions[0])
    Region: CirclePixelRegion
    center: PixCoord(x=2896.5, y=5056.5)
    radius: 381.9716

Serialisation is done using the `~regions.fits_region_objects_to_table`
function::

    >>> from regions import fits_region_objects_to_table
    >>> table_ouput = fits_region_objects_to_table(regions)
    >>> print(table_ouput)
         X [4]            Y [4]        SHAPE       R [4]         ROTANG  COMPONENT
          pix              pix                      pix           deg
    ---------------- ---------------- ------- ---------------- --------- ---------
       2896.5 .. 0.0    5056.5 .. 0.0  circle  381.9716 .. 0.0       0.0         1
    5282.0541 .. 0.0 4854.5699 .. 0.0  ROTBOX 1303.4597 .. 0.0 28.395178         2
       2944.5 .. 0.0    3472.5 .. 0.0 ellipse     288.0 .. 0.0  337.4048         3
        341.0 .. 0.0     345.0 .. 0.0  ROTBOX      56.0 .. 0.0      65.0         4
        341.0 .. 0.0     345.0 .. 0.0 ANNULUS      56.0 .. 0.0       0.0         5
        341.0 .. 0.0     345.0 .. 0.0   point       0.0 .. 0.0       0.0         6
        341.0 .. 0.0     345.0 .. 0.0   point       0.0 .. 0.0       0.0         7
          1.0 .. 4.0       5.0 .. 8.0 polygon       0.0 .. 0.0       0.0         8
         10.0 .. 0.0       5.5 .. 0.0  ROTBOX      10.0 .. 0.0       0.0         9

The `~regions.write_fits` and `~regions.read_fits_region`
functions write as well as read from a file in addition to doing the
region serialisation and parsing:

.. doctest-skip::

    >>> from regions import CirclePixelRegion, PixCoord, write_fits
    >>> reg_pixel = CirclePixelRegion(PixCoord(1, 2), 5)
    >>> write_fits('regions_output.fits', regions=[reg_pixel], overwrite=True)
