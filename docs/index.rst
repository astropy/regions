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

Defining Regions
----------------

.. toctree::
   :maxdepth: 1

   shapes
   metadata

Using Regions
-------------

.. toctree::
   :maxdepth: 1

   contains
   compound
   masks
   spherical_frame_transform
   spherical_bounding_regions
   plotting

I/O & Reference
---------------------

.. toctree::
   :maxdepth: 1

   region_io
   shapely
   api


Examples
========

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

The next example demonstrates how to plot spherical regions (e.g., a circle
and a longitude/latitude "range") on a full sky image (seen in the
`visualization demo here`_) using Matplotlib.

.. plot::

    from astropy.wcs import WCS
    from astropy.io import fits
    from astropy.utils.data import get_pkg_data_filename
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.visualization.wcsaxes.frame import EllipticalFrame
    from matplotlib import patheffects
    from matplotlib import pyplot as plt

    from regions import CircleSphericalSkyRegion, RangeSphericalSkyRegion

    image_file = get_pkg_data_filename('allsky/allsky_rosat.fits')
    image_data = fits.getdata(image_file, ext=0, memmap=False)
    image_hdr = fits.getheader(image_file)
    wcs = WCS(image_hdr)


    fig, ax = plt.subplots(figsize=(8,8),
                           subplot_kw=dict(projection=wcs,
                                           frame_class=EllipticalFrame))

    path_effects=[patheffects.withStroke(linewidth=3, foreground='black')]
    ax.coords.grid(color='white')
    ax.coords['glon'].set_ticklabel(color='white', path_effects=path_effects)
    ax.set_axisbelow(True)
    ax.set_xlabel("Galactic $\ell$", fontsize=14, labelpad=8)
    ax.set_ylabel("Galactic $b$", fontsize=14)
    im = ax.imshow(image_data, vmin=0., vmax=300., cmap='gray')

    # Clip the image to the frame
    im.set_clip_path(ax.coords.frame.patch)

    # Image is in galactic coordinates
    circ = CircleSphericalSkyRegion(SkyCoord(80,45,unit=u.deg,frame="galactic"),
                                    25*u.deg)
    circ.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":1000}
    ).plot(ax=ax, edgecolor='red',
           facecolor='none', lw=2)

    range = RangeSphericalSkyRegion(
        frame="galactic",
        longitude_range=[210,330]*u.deg,
        latitude_range=[-55,-15]*u.deg,
    )
    range.to_pixel(
       wcs=wcs,
       include_boundary_distortions=True,
       discretize_kwargs={"n_points":250}
    ).plot(ax=ax, edgecolor='blue',
           facecolor='none', lw=2)
