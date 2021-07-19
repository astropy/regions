Combining Regions
=================

There are a few ways to combine any two `~regions.Region` objects
into a compound region, i.e., a `~regions.CompoundPixelRegion` or
`~regions.CompoundSkyRegion` object.

Let's start by defining two sky regions::

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion

    >>> circle1 = CircleSkyRegion(
    ...    center=SkyCoord(1, 2, unit='deg', frame='galactic'),
    ...    radius=Angle('5 deg'))
    >>> circle2 = CircleSkyRegion(
    ...    center=SkyCoord(-4, 3, unit='deg', frame='galactic'),
    ...    radius=Angle('3 deg'))


Intersection (AND)
------------------

To create an intersection compound region, use either the ``&`` operator
or the :meth:`~regions.Region.intersection` method::

    >>> comp_region = circle1 & circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function and_>

    >>> comp_region = circle1.intersection(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function and_>


Union (OR)
----------

To create a union compound region, use either the ``|`` operator or the
:meth:`~regions.Region.union` method::

    >>> comp_region = circle1 | circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function or_>

    >>> comp_region = circle1.union(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function or_>


Symmetric Difference (XOR)
--------------------------

To create a symmetric difference compound region, use either the ``^``
operator or the :meth:`~regions.Region.symmetric_difference` method::

    >>> comp_region = circle1 ^ circle2
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function xor>

    >>> comp_region = circle1.symmetric_difference(circle2)
    >>> print(comp_region)
    Region: CompoundSkyRegion
    region1: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (1., 2.)>
    radius: 5.0 deg
    region2: Region: CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
        (356., 3.)>
    radius: 3.0 deg
    operator: <built-in function xor>


Example Illustrating Compound Regions
-------------------------------------

.. plot:: plot_compound.py
    :include-source: false
