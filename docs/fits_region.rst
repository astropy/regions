.. _gs-fits:

Reading/writing to FITS region files
====================================

The regions package provides the functionality to serialise and de-serialise
Python lists of `~regions.Region` objects to FITS region file. A `FITS
Region Binary Table <https://fits.gsfc.nasa.gov/registry/region.html>`_
defines a spatial region of a two-dimensional image in pixels.
The name of the header must be ``REGION`` for the `~regions.read_fits_region`
to parse the table. The file is read using the `~astropy.table.Table.read` method
which converts it into a `~astropy.table.Table` object.
`~regions.FITSRegionParser` parses this table. It converts the FITS region table to a
`~regions.ShapeList` object, which is a list of `~regions.Shape`, each
representing one FITS region in pixels. It should be noted that the regions are
actually described in sky coordinates. It is converted explicitly to sky regions
in the `~regions.read_fits_region` using the `to_sky` method.
The `~astropy.wcs.WCS` object is created with the help of FITS
header of the file. The `~regions.Shape` objects can be converted to
`~regions.Region` objects. Some of the regions such as ``rectangle``,
``rotrectangle``, ``pie``, ``sector`` are not supported by this
package.

.. code-block:: python

    >>> from regions import FITSRegionParser
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

    >>> parser = FITSRegionParser(table)
    >>> print(parser.shapes[0])
    Shape
    Type : reg
    Coord sys : physical
    Region type : circle
    Meta: {'tag': '1'}
    Composite: False
    Include: False

    >>> regions = parser.shapes.to_regions() # pixel regions
    >>> print(regions[0])
    Region: CirclePixelRegion
    center: PixCoord(x=2896.5, y=5056.5)
    radius: 381.9716

Serialisation is done using the `~regions.fits_region_objects_to_table` function

.. code-block:: python

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

The `~regions.write_fits_region` and `~regions.read_fits_region` functions
write as well as read from a file in addition to doing the region serialisation
and parsing.

.. code-block:: python

    >>> from regions import CirclePixelRegion, PixCoord, write_fits_region
    >>> reg_pixel = CirclePixelRegion(PixCoord(1, 2), 5)
    >>> write_fits_region('regions_output.fits', regions=[reg_pixel], overwrite=True)
