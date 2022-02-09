# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

from astropy.utils.exceptions import AstropyUserWarning

from .core import ds9_valid_symbols, DS9ParserError

__all__ = []


def _split_raw_metadata(raw_metadata):
    """
    Split the raw metadata into meta and visual dictionaries.

    The raw metadata is not changed, except unsupported metadata is
    ignored.

    Parameters
    ----------
    raw_metadata : dict
        The raw metadata as read from the DS9 region file/str.

    Returns
    -------
    meta, visual : dict
        The raw metadata and visual metadata dictionaries.
    """
    # not currently supported by regions
    unsupported_meta = ('line', 'vector', 'ruler', 'compass')

    # TODO: include text in visual keys to support text annotations for
    # all DS9 regions?
    ds9_visual_keys = ('color', 'dash', 'dashlist', 'fill', 'font', 'point',
                       'textangle', 'textrotate', 'width')

    meta = {}
    visual = {}
    for key, value in raw_metadata.items():
        if key in unsupported_meta:
            # ignore this special case because line=0 0 (no arrows) works
            if key == 'line' and '1' not in value:
                continue
            warnings.warn(f'DS9 meta "{key}={value}" is unsupported and '
                          'will be ignored', AstropyUserWarning)

        if key in ds9_visual_keys:
            visual[key] = value
        else:
            meta[key] = value

    # set the default plotting style to 'ds9' when parsing ds9 data
    visual['default_style'] = 'ds9'

    return meta, visual


def _translate_visual_metadata(shape, visual_meta):
    """
    Translate ds9 visual metadata to dictionary of matplotlib keywords.

    Parameters
    ----------
    shape : str
        The DS9 region shape (e.g., 'circle').

    visual_meta : dict
        The visual metadata.

    Returns
    -------
    visual : dict
        The visual dictionaries of region metadata translated
        to matplotlib keywords.
    """
    meta = visual_meta.copy()

    fill = meta.pop('fill', 0)
    # fill=1 is supported in DS9 only for the circle, ellipse, and box
    # shapes
    if fill == 1 and shape in ('circle', 'ellipse', 'box'):
        meta['fill'] = True

    dash = meta.pop('dash', 0)
    dashlist = meta.pop('dashlist', None)
    if int(dash) == 1:
        if dashlist is not None:
            dashes = tuple(int(i) for i in dashlist.split())
            meta['linestyle'] = (0, dashes)
        else:
            meta['linestyle'] = 'dashed'

        if shape == 'point':
            warnings.warn('dashed lines are unsupported for DS9 point '
                          'regions and will be ignored', AstropyUserWarning)
            meta.pop('linestyle', None)

    # "point=symbol [size]"; size is optional, symbol is not
    point = meta.pop('point', None)
    if point is not None:
        point_ = point.split()
        if len(point_) == 1:
            ds9_marker = point_[0]
        elif len(point_) == 2:
            ds9_marker, meta['markersize'] = point_
        else:
            raise ValueError(f'invalid point data "{point}"')
        meta['marker'] = ds9_valid_symbols[ds9_marker]

    font = meta.pop('font', None)
    if font is not None:
        (meta['fontname'], meta['fontsize'], meta['fontweight'],
         meta['fontstyle']) = font.split()

        # fontsize is a str
        try:
            meta['fontsize'] = int(meta['fontsize'])
        except ValueError:
            raise DS9ParserError('font size must be an integer, got '
                                 f'{meta["fontsize"]}') from None
        meta['fontstyle'] = meta['fontstyle'].replace('roman', 'normal')

    if shape == 'point':
        width = meta.pop('width', None)
        if width is not None:
            meta['markeredgewidth'] = width
    else:
        invalid = ('marker', 'markersize')
        for key in invalid:
            meta.pop(key, None)

    if shape == 'text':
        textangle = meta.pop('textangle', None)
        if textangle is not None:
            textrotate = meta.pop('textrotate', None)
            if textrotate != 0:  # rotate if None or 1
                meta['rotation'] = textangle

        # remove invalid mpl kwargs for matplotlib.text.Text
        invalid = ('linestyle', 'linewidth', 'fill')
        for key in invalid:
            meta.pop(key, None)
    else:
        # these kwarg are valid only for matplotlib.text.Text
        invalid = ('textangle', 'textrotate')
        for key in invalid:
            meta.pop(key, None)

    if shape not in ('point', 'line', 'text'):
        color = meta.pop('color', None)
        if color is not None:
            meta['facecolor'] = color
            meta['edgecolor'] = color

    return meta
