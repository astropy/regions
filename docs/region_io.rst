.. _regions_io:

Reading/Writing Region Files
****************************

The regions package provides a unified interface for reading, writing,
parsing, and serializing regions data in different standard formats.


Regions I/O
===========

The `~regions.Regions` class (which represents a list
of `~regions.Region` objects) includes four methods,
:meth:`~regions.Regions.read`, :meth:`~regions.Regions.write`,
:meth:`~regions.Regions.parse`, and :meth:`~regions.Regions.serialize`,
that make it possible to read, write, parse, and serialize region files
or region data.

The :class:`~regions.Regions` class has built-in support for various
input and output formats. A full list of the supported formats and is
shown in the table below. The ``Suffix`` column indicates the standard
filename suffix for a particular format.

============  ==============  ========================
   Format       Suffix           Description
============  ==============  ========================
       crtf   .crtf           `CASA Region Text Format <https://casadocs.readthedocs.io/en/stable/notebooks/image_analysis.html#Region-File-Format>`_
        ds9   .reg, .ds9      `DS9 Region Format <http://ds9.si.edu/doc/ref/region.html>`_
       fits   .fits           `FITS Region Binary Table <https://fits.gsfc.nasa.gov/registry/region.html>`_
============  ==============  ========================

Use the :meth:`~regions.Regions.get_formats` method to get the
registered I/O formats as a :class:`~astropy.table.Table`::

    >>> from regions import Regions
    >>> print(Regions.get_formats())
    Format Parse Serialize Read Write Auto-identify
    ------ ----- --------- ---- ----- -------------
      crtf   Yes       Yes  Yes   Yes           Yes
       ds9   Yes       Yes  Yes   Yes           Yes
      fits   Yes       Yes  Yes   Yes           Yes


Read
----

To read in a region file, first import the :class:`~regions.Regions`
class, then call the :meth:`~regions.Regions.read` method with the name
of the file and the file format, e.g.:

.. doctest-skip::

    >>> from regions import Regions
    >>> regions = Regions.read('my_regions.reg', format='ds9')

It is also possible to load regions directly from the internet using
URLs, e.g.:

.. doctest-skip::

    >>> regions = Regions.read('https://some.domain.edu/my_regions.reg', format='ds9')

If the ``format`` keyword is not set, the :meth:`~regions.Regions.read`
method will attempt to automatically determine the file format from the
filename extension or the file contents:

.. doctest-skip::

    >>> regions = Regions.read('my_regions.reg')

The underlying file handler will also automatically
detect various compressed data formats and transparently
uncompress them if supported by the Python installation (see
`~astropy.utils.data.get_readable_fileobj`):

.. doctest-skip::

    >>> regions = Regions.read('my_regions.reg.gz')


Write
-----

Use the :meth:`~regions.Regions.write` method to write regions to
a region file. Like the :meth:`~regions.Regions.read` method, the
filename extension will be used to automatically define the format if
unspecified.

.. doctest-skip::

    >>> regions.write('my_regions.crtf', format='crtf')
    >>> regions.write('my_regions.reg')

If the file already exists and you want to overwrite it, then set the
``overwrite`` keyword to `True`:

.. doctest-skip::

    >>> regions.write('my_regions.reg', overwrite=True)


Parse
-----

Region data in the form of a string or table may also be
parsed into a :class:`~regions.Regions` object by using the
:meth:`~regions.Regions.parse` method. The ``format`` keyword must be
specified when parsing region data. A region string can be parsed for
the ``crtf`` and ``ds9`` formats, while a `~astropy.table.Table` can be
parsed for the ``fits`` format::

    >>> from regions import Regions
    >>> regions_str = '# Region file format: DS9\nimage\ncircle(147.10,254.17,3.1) # color=green\n'
    >>> regions = Regions.parse(regions_str, format='ds9')
    >>> print(regions)
    [<CirclePixelRegion(center=PixCoord(x=146.1, y=253.17), radius=3.1)>]
    >>> print(regions[0].visual)
    {'default_style': 'ds9', 'facecolor': 'green', 'edgecolor': 'green'}


Serialize
---------

Regions can be serialized to a string or table by using the
:meth:`~regions.Regions.serialize` method. The ``format`` keyword must
be specified when serializing region data. Serializing regions to the
``crtf`` or ``ds9`` format will produce a region string, while the
``fits`` format will produce a `~astropy.table.Table`:

.. doctest-skip::

    >>> regions_str = regions.serialize(format='ds9')


Region Classes I/O
==================

Additionally, all of the Region classes (i.e., :class:`~regions.Region`
subclasses) support the :meth:`~regions.Regions.write` and
:meth:`~regions.Regions.serialize` methods.

As an example, let's create a :class:`~regions.CircleSkyRegion` object::

    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> from regions import CircleSkyRegion
    >>> coord = SkyCoord(202.469575, 47.19525833, unit='deg', frame='fk5')
    >>> region = CircleSkyRegion(coord, radius=0.01 * u.deg)
    >>> region
    <CircleSkyRegion(center=<SkyCoord (FK5: equinox=J2000.000): (ra, dec) in deg
    (202.469575, 47.19525833)>, radius=0.01 deg)>

To serialize the region::

    >>> region.serialize(format='ds9', precision=8)
    '# Region file format: DS9 astropy/regions\nj2000\ncircle(202.46957500,47.19525833,0.01000000)\n'

To write the region to a region file:

.. doctest-skip::

    >>> region.write('my_region.ds9', format='ds9')

Use the :meth:`~regions.Region.get_formats` method to list all available
formats and methods for the :class:`~regions.Region` subclasses::

    >>> print(region.get_formats())
    Format Parse Serialize Read Write Auto-identify
    ------ ----- --------- ---- ----- -------------
      crtf    No       Yes   No   Yes           Yes
       ds9    No       Yes   No   Yes           Yes
      fits    No       Yes   No   Yes           Yes


Region File Format Limitations
==============================

.. toctree::
    :maxdepth: 1

    ds9_io
    fits_io
