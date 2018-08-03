.. _gs-fits:

Reading/writing to FITS region files
====================================

The regions package provides the functionality to serialise and de-serialise
Python lists of `~regions.Region` objects to FITS region file. It is a FITS
Region Binary Table that defines a spatial region of a two-dimensional image in
pixels. The name of the header must be ``REGION`` for the
`~regions.read_fits_region` to parse the table. The file is read using the
`~astropy.table.Table.read` method which gets converted into an
`~astropy.table.Table` object. This table is then parsed using
`~regions.FITSRegionParser`. It converts the FITS region table to
`~regions.ShapeList` object, which is a list of `~regions.Shape` each
representing one FITS region. The `~regions.Shape` objects can be converted to
`~regions.Region` objects. Some of the regions such as ``rectangle``,
``rotrectangle``, ``pie``, ``sector`` are not supported by this
package.

.. code-block:: python

    >>> from regions import FITSRegionParser
    >>> from astropy.table import Table
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> filename = get_pkg_data_filename('data/regions.fits',
    ...                                  package='regions.io.fits.tests')
    >>> table = Table.read(filename)
    >>> print(table)
        X         Y      SHAPE          R [2]            ROTANG  COMPONENT
      pix       pix                     pix              deg
    --------- --------- ------- ---------------------- --------- ---------
       2896.5    5056.5  CIRCLE        381.9716 .. 0.0       0.0         1
    5282.0541 4854.5699     BOX 1303.4597 .. 655.09466 28.395178         2
       2944.5    3472.5 ELLIPSE         288.0 .. 512.0  337.4048         3
        341.0     345.0  ROTBOX           56.0 .. 78.0      65.0         4
        341.0     345.0 ANNULUS           56.0 .. 78.0       0.0         5
        341.0     345.0   POINT             0.0 .. 0.0       0.0         6
        341.0     345.0                     0.0 .. 0.0       0.0         7
    >>> parser = FITSRegionParser(table)
    >>> print(parser.shapes[0])
    Shape
    Type : reg
    Coord sys : physical
    Region type : circle
    Meta: {'tag': '1'}
    Composite: False
    Include: False
    >>> regions = parser.shapes.to_regions()
    >>> print(regions[0])
    Region: CirclePixelRegion
    center: PixCoord(x=2896.5, y=5056.5)
    radius: 381.9716
    >>> print(regions[0].meta)
    {'tag': '1', 'include': False}

Serialisation is done using the `~regions.fits_region_objects_to_table` function

.. code-block:: python

    >>> from regions import fits_region_objects_to_table
    >>> table_ouput = fits_region_objects_to_table(regions)
    >>> print(table_ouput)
      X [1]     Y [1]    SHAPE       R [4]         ROTANG  COMPONENT
       pix       pix                  pix           deg
    --------- --------- ------- ---------------- --------- ---------
       2896.5    5056.5  circle  381.9716 .. 0.0       0.0         1
    5282.0541 4854.5699  ROTBOX 1303.4597 .. 0.0 28.395178         2
       2944.5    3472.5 ellipse     288.0 .. 0.0  337.4048         3
        341.0     345.0  ROTBOX      56.0 .. 0.0      65.0         4
        341.0     345.0 ANNULUS      56.0 .. 0.0       0.0         5
        341.0     345.0   point       0.0 .. 0.0       0.0         6
        341.0     345.0   point       0.0 .. 0.0       0.0         7


There's also `~regions.write_fits_region` and `~regions.read_fits_region` to
write as well as read from a file in addition to doing the region serialisation
and parsing.

.. code-block:: python

    >>> from astropy.io import fits
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> from regions import read_fits_region, write_fits_region
    >>> file_read = get_pkg_data_filename('data/region.fits',
    ...                                   package='regions.io.fits.tests')
    >>> hdul = fits.open(file_read)
    >>> regions = read_fits_region(file_read)
    >>> filename = 'region_ouput.fits'
    >>> write_fits_region(regions, hdul[1].header, filename)
