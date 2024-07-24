.. include:: references.txt

###############
Astropy Regions
###############

**regions** is an in-development `coordinated package`_ of `Astropy`_
for region handling.

To get an overview of available features, see :ref:`getting_started`.

The eventual goal is to merge this package into the Astropy core
package.

* Code : `GitHub repository`_
* Contributors : https://github.com/astropy/regions/graphs/contributors


Getting Started
===============

.. toctree::
   :maxdepth: 1

   install
   getting_started
   contributing
   license
   changelog


User Documentation
==================

.. toctree::
   :maxdepth: 1

   shapes
   metadata
   contains
   compound
   masks
   plotting
   region_io
   shapely
   api


Example
=======

The following example shows how to read a DS9 region file and plot the
regions on an image using Matplotlib.

.. plot::

    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    from matplotlib import pyplot as plt

    from regions import Regions

    image_file = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
    image_data = fits.getdata(image_file, ext=0, memmap=False)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(image_data, cmap='gray')
    #ax.set_ylim([-0.5, 892.5])

    region_file = get_pkg_data_filename('data/plot_image.reg',
                                        package='regions.io.ds9.tests')
    regions = Regions.read(region_file, format='ds9')
    for i, region in enumerate(regions):
        region.plot(ax=ax)
