Checking for Points Inside Regions
==================================

Points Inside Planar Regions
----------------------------

Let's start by defining both a planar sky and pixel region::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion, PixCoord, CirclePixelRegion

    >>> sky_center = SkyCoord(42, 43, unit='deg')
    >>> sky_radius = Angle(25, 'deg')
    >>> sky_region = CircleSkyRegion(sky_center, sky_radius)
    >>> print(sky_region)
    Region: CircleSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 25.0 deg

    >>> pixel_center = PixCoord(x=42, y=43)
    >>> pixel_radius = 42
    >>> pixel_region = CirclePixelRegion(pixel_center, pixel_radius)
    >>> print(pixel_region)
    Region: CirclePixelRegion
    center: PixCoord(x=42, y=43)
    radius: 42

Let's also define a WCS object using our example dataset::

    >>> from regions import make_example_dataset
    >>> dataset = make_example_dataset(data='simulated')
    >>> wcs = dataset.wcs

To test if a given point is inside or outside the region, the Python
``in`` operator can be called::

    >>> from regions import PixCoord
    >>> PixCoord(55, 40) in pixel_region
    True
    >>> PixCoord(55, 200) in pixel_region
    False

The ``in`` operator works only for scalar coordinates and pixel regions.
If you try to use ``in`` for non-scalar coordinates, you'll get a
`ValueError`::

    >>> pixcoord = PixCoord([50, 50], [10, 60])
    >>> pixcoord in pixel_region
    Traceback (most recent call last):
    ...
    ValueError: coord must be scalar. coord=PixCoord(x=[50 50], y=[10 60])

If you have arrays of coordinates, use the `regions.SkyRegion.contains`
or `regions.PixelRegion.contains` methods::

    >>> pixcoords = PixCoord.from_sky(sky_center, wcs)
    >>> pixel_region.contains(pixcoords)
    True

Note that `regions.SkyRegion.contains` requires a WCS to be passed::

    >>> skycoord = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sky_region.contains(skycoord, wcs)
    array([False, True])


Points Inside Spherical Regions
-------------------------------

For `~regions.SphericalSkyRegion` objects, checking whether point(s) are
contained inside that region requires no other input --- since these
regions are defined with a the spherical geometry, and not a projected geometry
(as captured through the projection encoded in a WCS) as in
`~regions.SkyRegion`.

Let's define a spherical sky region::

    >>> from regions import CircleSphericalSkyRegion

    >>> sph_sky_center = SkyCoord(42, 43, unit='deg')
    >>> sph_sky_radius = Angle(25, 'deg')
    >>> sph_sky_region = CircleSphericalSkyRegion(sph_sky_center,
    ...                                           sph_sky_radius)
    >>> print(sph_sky_region)
    Region: CircleSphericalSkyRegion
    center: <SkyCoord (ICRS): (ra, dec) in deg
        (42., 43.)>
    radius: 25.0 deg

Use the `~regions.SphericalSkyRegion.contains()` method to determine which
point(s) lie inside or outside the region::

    >>> skycoord = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sph_sky_region.contains(skycoord)
    array([False, True])
