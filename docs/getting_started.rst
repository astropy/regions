===============
Getting started
===============

Let's make a region:

.. code-block:: python

    from astropy.coordinates import Angle, SkyCoord
    import regions

    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = regions.shapes.CircleSkyRegion(center, radius)
    print(region)


TODO: do more things ...
