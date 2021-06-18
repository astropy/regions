.. _unified-io:

Reading/Writing Region Files
****************************

The regions package provides a unified interface for reading, writing,
parsing, and serializing regions data in different formats.


Getting Started with Regions I/O
================================

The `~regions.Regions` class includes four methods,
:meth:`~regions.Regions.read`, :meth:`~regions.Regions.write`,
:meth:`~regions.Regions.parse`, and :meth:`~regions.Regions.serialize`,
that make is possible to read, write, parse, and serialize region files
or region data.

The :class:`~regions.Regions` class has built-in support for various
input and output formats. A full list of the supported formats and is
shown in the table below. The ``Suffix`` column indicates the standard
filename suffix for a particular format.

============  ==============  ========================
   Format       Suffix           Description
============  ==============  ========================
       crtf   .crtf           `CASA Region Text Format <https://casa.nrao.edu/casadocs/casa-6.1.0/imaging/image-analysis/region-file-format>`_
        ds9   .reg, .ds9      `DS9 Region Format <http://ds9.si.edu/doc/ref/region.html>`_
       fits   .fits           `FITS Region Binary Table <https://fits.gsfc.nasa.gov/registry/region.html>`_
============  ==============  ========================


Examples
--------

Read
^^^^

To read in a region file, first import the :class:`~regions.Regions`
class, then call the :meth:`~regions.Regions.read` method with the name
of the file and the file format, e.g.:

.. doctest-skip::

    >>> from regions import Regions
    >>> regions = Regions.read('my_regions.reg', format='ds9')

It is also possible to load tables directly from the internet using
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
^^^^^

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
^^^^^

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
    {'color': 'green'}


Serialize
^^^^^^^^^

Regions can be serialized to a string or table by using the
:meth:`~regions.Regions.serialize` method. The ``format`` keyword must
be specified when serializing region data. Serializing regions to the
``crtf`` or ``ds9`` format will produce a region string, while the
``fits`` format will produce a `~astropy.table.Table`:

.. doctest-skip::

    >>> regions_str = regions.serialize(format='ds9')
