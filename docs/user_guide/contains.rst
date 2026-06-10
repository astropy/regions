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
    np.True_

Note that `regions.SkyRegion.contains` requires a WCS to be passed::

    >>> skycoord = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sky_region.contains(skycoord, wcs)
    array([False, True])


Points Inside Spherical Regions
-------------------------------

For `~regions.SphericalSkyRegion` objects, checking whether point(s) are
contained inside that region requires no other input --- since these
regions are defined with a spherical geometry, and not a projected geometry
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

Use the `regions.SphericalSkyRegion.contains` method to determine which
point(s) lie inside or outside the region::

    >>> skycoord = SkyCoord([50, 50], [10, 60], unit='deg')
    >>> sph_sky_region.contains(skycoord)
    array([False, True])


Boundary Behavior: ``contains`` versus ``covers``
-------------------------------------------------

The ``contains`` method excludes the region boundary: a point that lies
exactly on an edge or vertex is considered to be *outside* the region.
This is consistent with `Shapely's
<https://shapely.readthedocs.io/>`_ ``contains`` predicate and the
`DE-9IM <https://en.wikipedia.org/wiki/DE-9IM>`_ semantics used by most
GIS tools.

If you instead want boundary points to be treated as *inside* the
region, use the ``covers`` method. It is available on pixel regions
(`regions.PixelRegion.covers`), planar sky regions
(`regions.SkyRegion.covers`), and spherical sky regions
(`regions.SphericalSkyRegion.covers`), and takes the same arguments as
the corresponding ``contains`` method.

For example, the point ``(84, 43)`` lies exactly on the boundary of our
circular pixel region (centered at ``(42, 43)`` with a radius of
``42``), so it is excluded by ``contains`` but included by ``covers``::

    >>> boundary_point = PixCoord(84, 43)
    >>> pixel_region.contains(boundary_point)
    np.False_
    >>> pixel_region.covers(boundary_point)
    np.True_

For points strictly inside or strictly outside the region, the two
methods return identical results::

    >>> pixel_region.contains(PixCoord(42, 43))
    np.True_
    >>> pixel_region.covers(PixCoord(42, 43))
    np.True_

.. note::

    The ``covers`` method is currently implemented for circle, ellipse,
    rectangle, polygon, and annulus regions (and their sky and
    spherical-sky equivalents); calling it on a region type that does
    not support it raises a `NotImplementedError`.
