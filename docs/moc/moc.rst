.. _moc:

Multi-Order Coverage maps
=========================

.. _moc-intro:

Introduction
------------

A MOC describes the coverage of an arbitrary region on the unit sphere.
MOCs are usually used for describing the global coverage of catalog/image surveys such as GALEX or SDSS.
A MOC corresponds to a list of `HEALPix <https://healpix.sourceforge.io/>`__ cells at different depths.
This class gives you the possibility to:

1. Define `~regions.MOCSkyRegion` objects:

- From a FITS file that stores HEALPix cells (see `~regions.MOCSkyRegion.from_fits`). :ref:`[1] <from_fits>`
- Directly from a list of HEALPix cells expressed either as a numpy structural array (see `~regions.MOCSkyRegion.from_cells`) or a simple
    python dictionnary (see `~regions.MOCSkyRegion.from_json`).
- From a list of sky coordinates (see `~regions.MOCSkyRegion.from_skycoord`, `~regions.MOCSkyRegion.from_lonlat`).
- From a convex/concave polygon (see `~regions.MOCSkyRegion.from_polygon`). :ref:`[3] <from_polygon>`
- From a cone (will be implemented in a next version).

2. Perform fast logical operations between `~regions.MOCSkyRegion` objects: :ref:`[4] <logical_operations>`

- The `~regions.MOCSkyRegion.intersection` 
- The `~regions.MOCSkyRegion.union`
- The `~regions.MOCSkyRegion.difference`
- The `~regions.MOCSkyRegion.complement`

3. :ref:`plot_moc`:

- Draw the MOC with its HEALPix cells (see `~regions.MOCSkyRegion.fill`)
- Draw the perimeter of a MOC (see `~regions.MOCSkyRegion.border`)

4. Get the sky coordinates defining the border(s) of `~regions.MOCSkyRegion` objects (see `~regions.MOCSkyRegion.get_boundaries`).

5. Serialize `~regions.MOCSkyRegion` objects to `astropy.io.fits.HDUList` or JSON dictionary and save it to a file. :ref:`[6] <serialize>`

More informations about MOCs can be found by reading the following `IVOA paper <http://www.ivoa.net/documents/MOC/20140602/REC-MOC-1.0-20140602.pdf>`_.

Examples
--------

.. _from_fits:

Load a `~regions.MOCSkyRegion` from a FITS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A MOC is often loaded from a FITS file. A FITS file describing a MOC contains an `astropy.io.fits.BinTableHDU`
table where one can find a list of 64 bits unsigned numbers. Each of these numbers (usually called uniq numbers)
describes one HEALPix cell at a specific depth.

For the purpose of the example, we will load the MOC describing the coverage
of the ``P-GALEXGR6-AIS-FUV`` sky survey. It is a rather complex MOC containing lots of holes
and spread over a good range of different depths.

.. code-block:: python

    >>> from regions import MOCSkyRegion
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))

It is also possible to create a MOC from a list of `~astropy.coordinates.SkyCoord`, lon and lat `~astropy.units.Quantity`,
a python dict storing HEALPix cells etc... Please refer to the `~regions.MOCSkyRegion` API for more details.


.. _plot_moc:

Plot a `~regions.MOCSkyRegion`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example:

- Load the coverage of GALEX from a FITS file.
- Plot the MOC by:

1. Defining a matplotlib figure.
2. Defining an astropy WCS representing the field of view of the plot.
3. Create a `~regions.MOCPixelRegion` from the moc sky region instance and plot it.
4. Set the axis labels, a title, enable the grid and plot the final figure.

.. plot:: moc/plot_moc.py
    :include-source: true

It is also possible to use the `regions.MOCSkyRegion.fill` and `regions.MOCSkyRegion.border` methods from `MOCSkyRegion` without casting it to a `MOCPixelRegion`.
These methods take an `~astropy.wcs.WCS` and a `matplotlib.axes.Axes` to draw the MOC into it.

.. plot:: moc/plot_moc_no_cast.py
    :include-source: true

.. _from_polygon:

Create `~regions.MOCSkyRegion` from a polygon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example creates a `~regions.MOCSkyRegion` instance from a list of sky coordinates.
These sky coordinates refer to the vertices of a polygon.
Concave/Convex polygons are accepted but not self intersecting ones.
Two methods allow you to create a `~regions.MOCSkyRegion` instance from a polygon:

- `~regions.MOCSkyRegion.from_polygon` asks for two `~astropy.units.Quantity` storing the longitudes (resp. the latitudes) of the polygon vertices.
- `~regions.MOCSkyRegion.from_polygon_skycoord` asks for an `~astropy.coordinates.SkyCoord` storing the vertices of the polygon.

.. plot:: moc/plot_polygon.py
    :include-source: true

.. _logical_operations:

Intersection between GALEX and SDSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example:

- Load the coverages of SDSS and GALEX from FITS files.

.. code-block:: python

    >>> from regions import MOCSkyRegion
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))
    >>> sdss9 = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-SDSS9-r.fits', package='regions'))

- Compute their intersection

.. code-block:: python

    >>> intersection = galex.intersection(sdss9)

- Compute their union

.. code-block:: python

    >>> union = galex.union(sdss9)

- Plot the resulting intersection and union on a same matplotlib axis.

.. plot:: moc/plot_logical_operations.py
    :include-source: false

Check whether sky coordiantes fall inside a `~regions.MOCSkyRegion`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given ra and dec `~astropy.units.Quantity` instances, one can check if these sky positions
lie inside a `~regions.MOCSkyRegion`.

.. code-block:: python

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))
    >>> # Generate a thousand random positions on the sky 
    ... ra = np.random.randint(low=0, high=360, size=1000) * u.deg
    ... dec = np.random.randint(low=-90, high=90, size=1000) * u.deg
    >>> inside_mask = galex.contains(ra, dec)

.. plot:: moc/plot_contains.py
    :include-source: false

.. _serialize:

Write a `~regions.MOCSkyRegion` to a FITS/JSON file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can serialize a `~regions.MOCSkyRegion` into a `~astropy.io.fits.hdu.hdulist.HDUList` object or a
simple python dictionary (JSON format). The following block of code 
serializes the galex `~regions.MOCSkyRegion` instance defined in :ref:`[1] <from_fits>` to a FITS hdulist.

.. code-block:: python

    >>> hdulist = galex.serialize()
    >>> hdulist.info()
    Filename: (No file associated with this HDUList)
    No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU       4   ()      
    1                1 BinTableHDU     15   71002R x 1C   ['1J']   

The hdulist contains two tables, a primary and a binary one. The MOC is stored in the second one.
It consists of a list of the uniq numbers representing the set of HEALPix cells contained in the MOC (hence the 1d table).

If you want to write it, just call the `~regions.MOCSkyRegion.write` method with the path
you wish to save the MOC instance to.
