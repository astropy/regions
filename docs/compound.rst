
.. _gs-compound:

Combining regions
=================

There's a few ways to combine any two `~regions.Region` objects into a compound region,
i.e. a `~regions.CompoundPixelRegion` or `~regions.CompoundSkyRegion` object.

* The ``&`` operator calls the ``__and__`` method which calls the :meth:`~regions.Region.intersection` method
  to create an intersection compound region.
* The ``|`` operator calls the ``__or__`` method which calls the :meth:`~regions.Region.union` method
  to create a union compound region.
* The ``^`` operator calls the ``__xor__`` method which calls the :meth:`~regions.Region.symmetric_difference` method
  to create a symmetric difference compound region.

.. code-block:: python

    >>> from astropy.coordinates import Angle, SkyCoord
    >>> from regions import CircleSkyRegion
    >>> circle1 = CircleSkyRegion(
    ...    center=SkyCoord(1,2, unit='deg', frame='galactic'),
    ...    radius=Angle('5 deg')
    ...    )
    >>> circle2 = CircleSkyRegion(
    ...    center=SkyCoord(-4,3, unit='deg', frame='galactic'),
    ...    radius=Angle('3 deg'),
    ...    )
    >>> print(circle1 & circle2)
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
    >>> print(circle1 ^ circle2)
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

.. plot:: plot_compound.py
    :include-source: false
