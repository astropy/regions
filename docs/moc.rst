.. _moc:

Multi-Order Coverage maps
=========================

.. _moc-intro:

Introduction
------------

MOC stands for Multi-Order Coverage map. It is a new standard
for representing all-sky coverages. MOCs are defined as a list of HEALPix 
cells at different levels. Shallow HEALPix cells allow to describe big areas
of the sky while deeper ones are a good fit for representing smaller
regions. With that in mind, MOCs can describe the coverage of any sky survey.
The deepest level a MOC can be is limited to 29 which corresponds to a resolution of ~393.2uas.

The internal representation of a MOC uses the HEALPix nested numbering scheme for always keeping
the MOC consistent : 4 HEALPix neighbours belonging to the MOC (i.e. of indexes i, i+1, i+2, i+3 and of the same level) can be removed from it and replaced 
by the i//4's HEALPix cell of the level above (where // indicates the division operator followed by flooring the result).

More information about MOCs can be found by reading the following `IVOA paper <http://www.ivoa.net/documents/MOC/20140602/REC-MOC-1.0-20140602.pdf>`_.

Examples
--------

Loading from a FITS file
~~~~~~~~~~~~~~~~~~~~~~~~

A MOC is often loaded from a FITS file. A FITS file describing a MOC contains a `astropy.io.fits.BinTableHDU`
table where one can find a list of 64 bits numbers. Each of these numbers (usually called uniq numbers)
describes the index of one HEALPix cell at a specific level.

For the purpose of the example, we will load the MOC describing the coverage
of the ``P-GALEXGR6-AIS-FUV`` sky survey. It is a rather complex MOC containing lots of holes
and spread over a good range of different levels.

.. code-block:: python

    >>> from regions import MOCSkyRegion
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))

It is also possible to create a MOC from a list of `~astropy.coordinates.SkyCoord`, lon and lat `~astropy.units.Quantity`,
a python dict storing HEALPix cells etc... Please refer to the `~regions.MOCSkyRegion` API for more details.


Plot a `~regions.MOCSkyRegion`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's see how to plot it. For that purpose we need a `~astropy.wcs.WCS` instance representing the way 
our MOC will be projeted. We will initialize a simple WCS describing a galactical aitoff projection.

.. code-block:: python
    
    >>> import numpy as np
    >>> from astropy import wcs
    >>> w = wcs.WCS(naxis=2)
    ... w.wcs.crpix = [0, 0]
    ... w.wcs.cdelt = np.array([-5, 5])
    ... w.wcs.crval = [0, 0]
    ... w.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]


Once the WCS is defined, we perform the projection of the `~regions.MOCSkyRegion` according to the wcs object.
We get a `~regions.MOCPixelRegion` and next we can call its plot method to draw it 
using some matplotlib keyword arguments.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={'projection': w})
    >>> galex.to_pixel(w).plot(ax=ax, edgecolor='salmon', lw=1, facecolor='royalblue')
    >>> plt.title('P-GALEXGR6-AIS-FUV')
    ... plt.axis('equal')
    ... plt.xlabel('ra')
    ... plt.ylabel('dec')
    ... plt.xlim(-20, 20)
    ... plt.ylim(-10, 10)
    ... plt.grid(True)
    ... plt.show()


.. plot::

    from regions import MOCSkyRegion
    from astropy.utils.data import get_pkg_data_filename
    moc = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))

    # Configure a WCS
    import numpy as np
    from astropy import wcs
    w = wcs.WCS(naxis=2)

    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = np.array([-5, 5])
    w.wcs.crval = [0, 0]
    w.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={'projection': w})
    moc.to_pixel(w).plot(ax=ax, edgecolor='salmon', lw=1, facecolor='royalblue')

    plt.axis('equal')
    plt.title('P-GALEXGR6-AIS-FUV')
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.xlim(-30, 30)
    plt.ylim(-30, 30)
    plt.grid(True)
    plt.show()

Common operations between `~regions.MOCSkyRegion` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can perform the intersection, union and difference between two `~regions.MOCSkyRegion` objects.
Let is do the intersection between the coverage of ``P-GALEXGR6-AIS-FUV`` and the one of ``P-SDSS9-r``.

.. code-block:: python

    >>> from regions import MOCSkyRegion
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))
    >>> sdss9 = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-SDSS9-r.fits', package='regions'))

.. code-block:: python

    >>> result = galex.intersection(sdss9)

.. plot::

    from regions import MOCSkyRegion
    from astropy.utils.data import get_pkg_data_filename

    galex = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))
    sdss9 = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-SDSS9-r.fits', package='regions'))
    result = galex.intersection(sdss9)

    # Configure a WCS
    import numpy as np
    from astropy import wcs
    w = wcs.WCS(naxis=2)

    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = np.array([-5, 5])
    w.wcs.crval = [0, 0]
    w.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]

    def plot(ax, moc, title):
        moc.to_pixel(w).plot(ax=ax, edgecolor='salmon', lw=1, facecolor='royalblue')
        ax.axis('equal')
        ax.set_title(title)
        ax.set_xlabel('ra')
        ax.set_ylabel('dec')
        ax.set_xlim(-30, 30)
        ax.set_ylim(-30, 30)
        ax.grid(True)

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(30, 10), subplot_kw={'projection': w})

    plot(axes[0], galex, 'P-GALEXGR6-AIS-FUV')
    plot(axes[1], sdss9, 'P-SDSS9-r')
    plot(axes[2], result, 'P-GALEXGR6-AIS-FUV & P-SDSS9-r')

    plt.show()

Check whether sky positions fall inside a `~regions.MOCSkyRegion`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given ra and dec `astropy.units.Quantity` instances, one can check if these sky positions
lies inside a `~regions.MOCSkyRegion`.

.. code-block:: python

    >>> import numpy as np
    >>> from astropy import units as u
    >>> # Generate a thousand random positions on the sky 
    ... ra = np.random.randint(low=0, high=360, size=1000) * u.deg
    ... dec = np.random.randint(low=-90, high=90, size=1000) * u.deg
    >>> inside_mask = galex.contains(ra, dec)

.. plot::

    from regions import MOCSkyRegion
    from astropy.utils.data import get_pkg_data_filename
    moc = MOCSkyRegion.from_fits(get_pkg_data_filename('shapes/tests/data/P-GALEXGR6-AIS-FUV.fits', package='regions'))

    # Galactical aitoff projection
    import numpy as np
    from astropy import wcs
    w = wcs.WCS(naxis=2)

    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = np.array([-5, 5])
    w.wcs.crval = [0, 0]
    w.wcs.ctype = ["GLON-AIT", "GLAT-AIT"]

    # Create 1000 random positions on the sky
    from astropy import units as u
    ra = np.random.randint(low=0, high=360, size=1000) * u.deg
    dec = np.random.randint(low=-90, high=90, size=1000) * u.deg

    # Compute the mask of the positions inside the MOC
    inside_mask = moc.contains(ra, dec)

    # Project the positions with the WCS
    from astropy.wcs.utils import skycoord_to_pixel
    from astropy.coordinates import SkyCoord

    coords_inside = SkyCoord(ra[inside_mask], dec[inside_mask])
    coords_outside = SkyCoord(ra[~inside_mask], dec[~inside_mask])

    x_in, y_in = skycoord_to_pixel(coords_inside, wcs=w)
    x_out, y_out = skycoord_to_pixel(coords_outside, wcs=w)

    # Matplotlib code
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), subplot_kw={'projection': w})
    # Plot the MOC
    moc.to_pixel(w).plot(ax=ax, edgecolor='salmon', lw=1, facecolor='royalblue', zorder=1)

    plt.axis('equal')
    plt.title('Coordinates inside P-GALEXGR6-AIS-FUV')
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.grid(True)
    # Plot the projeted positions inside (resp. outside) the MOC
    plt.scatter(x_out, y_out, s=16, c='black', marker='^', zorder=2, label='outside')
    plt.scatter(x_in, y_in, s=16, c='red', marker='^', zorder=3, label='inside')
    plt.legend(loc='upper left')
    plt.show()

Write a `~regions.MOCSkyRegion` to a FITS/JSON file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can serialize a `~regions.MOCSkyRegion` to a `~astropy.io.fits.hdu.hdulist.HDUList` object or a
simple python dictionary (JSON format). The following block of code 
serializes the galex `~regions.MOCSkyRegion` object defined above to a FITS hdulist.

.. code-block:: python

    >>> hdulist = galex.serialize()
    >>> hdulist.info()
    Filename: (No file associated with this HDUList)
    No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU       4   ()      
    1                1 BinTableHDU     15   71002R x 1C   ['1J']   

The hdulist contains two tables, a primary and a binary one. The MOC is stored into the second one.
It consists of the list of uniq numbers representing the set of HEALPix cell contained in the MOC (hence the 1d table).

If you want to write it, just call the `~regions.MOCSkyRegion.write` method with the path
you wish to save the `~regions.MOCSkyRegion` instance to.

