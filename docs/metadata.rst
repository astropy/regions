.. include:: references.txt

Region Metadata
===============

A :class:`~regions.Region` object has both ``meta`` and
``visual`` attributes that store the metadata of the region. The
:class:`~regions.RegionMeta` and :class:`~regions.RegionVisual` classes
handle the general metadata and visual attributes, respectively. Both
behave like Python dictionaries.


General Metadata
----------------

The :class:`~regions.Region` ``meta`` attribute stores additional
metadata about the region such as labels, tags, comments, name, etc.
that are used for non-display tasks. It can also store the spectral
dimensions of the region.

The valid keys in the :class:`~regions.RegionMeta` object are:

* ``label``:

  - CRTF, DS9 (text label for a region)
  - Ex: meta['label'] = 'this is a circle'

* ``tag``:

  - DS9 (All regions may have zero or more tags associated with it,
    which may be used for grouping and searching.)

  - Ex: meta['tags'] = ['{Group 1}', '{Group 2}']}

* ``include``:

  - CRTF, DS9 (Region inclusion)
  - Possible Values: True, False
  - Ex: meta['include'] = True

* ``frame``:

  - CRTF (Frequency/Velocity Axis)
  - Possible values: 'REST', 'LSRK', 'LSRD', 'BARY', 'GEO', 'TOPO',
    'GALACTO', 'LGROUP', 'CMB'
  - Default: image value
  - Ex: meta['frame'] = 'TOPO'

* ``range``:

  - CRTF (Frequency/Velocity Range)
  - Possible units: GHz, MHz, kHz, km/s, Hz, channel, chan (=channel)
  - Default: image range
  - Format: [min, max]
  - Ex: meta['range'] = [-320 * u.m/u.sec, -330 * u.m/u.s]

* ``veltype``:

  - CRTF (Velocity Calculation)
  - Possible values: 'RADIO', 'OPTICAL', 'Z', 'BETA', 'GAMMA'
  - Default: image value
  - Ex: meta['veltype'] = 'RADIO'

* ``restfreq``:

  - CRTF (Rest Frequency)
  - Possible values: `~astropy.units.Quantity` object
  - Default: image value
  - Ex: meta['restfreq'] = Quantity("1.42GHz")

* ``corr``:

  - CRTF (Correlational Axis)
  - Possible values: 'I', 'Q', 'U', 'V', 'RR', 'RL', 'LR', 'LL',
    'XX', 'XY', 'YX', 'YY', 'RX', 'RY', 'LX', 'LY', 'XR', 'XL', 'YR',
    'YL', 'PP', 'PQ', 'QP', 'QQ', 'RCircular', 'LCircular', 'Linear',
    'Ptotal', 'Plinear', 'PFtotal', 'PFlinear', 'Pangle'
  - Default: all planes present in image
  - Ex: meta['corr'] = ['X', 'Y']

* ``comment``:

  - DS9, CRTF (Comment on the region)
  - Ex: meta['comment'] = 'Any comment for the region'

* ``line``:

  - DS9 (The line region may be rendered with arrows, one at each end.
    To indicate arrows, use the line property. A '1' indicates an arrow,
    '0' indicates no arrow.)
  - Ex: meta['line'] = [1, 1]

* ``name``

* ``select``

* ``highlite``

* ``fixed``

* ``edit``:

  - DS9 (The Edit property specifies if the user is allowed to edit the
    region via the GUI.)
  - Ex: meta['edit'] = 1

* ``move``:

  - DS9 (The Move property specifies if the user is allowed to move the
    region via the GUI. )
  - Ex: meta['move'] = 1

* ``rotate``:

  - DS9 (The Rotate property specifies if the user is allowed to rotate
    the region via the GUI. )
  - Ex: meta['rotate'] = 1

* ``delete``:

  - DS9 (The Delete property specifies if the user is allowed to delete
    the region via the GUI. )
  - Ex: meta['delete'] = 1

* ``source``

* ``background``


Visual Metadata
---------------

The :class:`~regions.Region` ``visual`` attribute stores visual
properties for the region such as color, font style, font size, line
width, line style, etc. The visual attributes are metadata meant to be
used to visualize regions, especially used by plotting libraries such as
`Matplotlib`_ .

The valid keys in the `~regions.RegionVisual` class are:

* ``color``: CRTF, DS9 (Region, symbol and text color)

  - Possible values: any color recognized by `Matplotlib`_, including
    hex values
  - Default: color=green
  - Ex: visual['color'] = 'blue'

* ``dash``: Render region using dashed lines using current dashlist
  value.

* ``font``: Name of the font.

* ``dashlist``: Sets dashed line parameters. This does not render the
  region in dashed lines.

* ``symsize``: Size of the symbol.

* ``symthick``: Thickness of the symbol.

* ``fontsize``: Size of the font.

* ``fontstyle``: Style of the font.

* ``usetex``: Boolean value whether the label uses TeX.

* ``labelpos``: Position of the label.

* ``labeloff``: Label offset.

* ``linewidth``: Width of the line.

* ``linestyle``: Style of the line.

* ``fill``: Boolean value whether the region is filled.

* ``line``: The line region may be rendered with arrows, one at each
  end. To indicate arrows, use the line property. A '1' indicates an
  arrow, '0' indicates no arrow.

* ``symbol``/``point``: CRTF, DS9 (Symbol for which a point region is
  described)

  - Ex: meta['symbol'] = 'point marker'
