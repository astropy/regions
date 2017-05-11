
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
    >>> type(circle1 & circle2)
    regions.core.compound.CompoundSkyRegion
    >>> print(circle1 ^ circle2)
    (CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
    (1.0, 2.0)>
    radius: 5.0 deg
    <built-in function xor>
    CircleSkyRegion
    center: <SkyCoord (Galactic): (l, b) in deg
    (356.0, 3.0)>
    radius: 3.0 deg)

.. plot:: plot_compound.py
    :include-source: false
