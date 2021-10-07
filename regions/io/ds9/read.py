# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import itertools
import re
import string
import warnings
from dataclasses import dataclass

from astropy.coordinates import Angle, frame_transform_graph
import astropy.units as u
from astropy.utils.data import get_readable_fileobj

from ...core import Regions, RegionMeta, RegionVisual
from ...core.registry import RegionsRegistry
from ..core import _Shape, _ShapeList, reg_mapping
from .core import (DS9RegionParserError, DS9RegionParserWarning,
                   valid_symbols_ds9)

from ...shapes import (CirclePixelRegion, CircleSkyRegion,
                       EllipsePixelRegion, EllipseSkyRegion,
                       RectanglePixelRegion, RectangleSkyRegion,
                       PolygonPixelRegion, PolygonSkyRegion,
                       RegularPolygonPixelRegion,
                       CircleAnnulusPixelRegion, CircleAnnulusSkyRegion,
                       EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion,
                       RectangleAnnulusPixelRegion, RectangleAnnulusSkyRegion,
                       LinePixelRegion, LineSkyRegion,
                       PointPixelRegion, PointSkyRegion,
                       TextPixelRegion, TextSkyRegion)


__all__ = []

# Regular expression to extract meta attributes
meta_regex = re.compile('([a-zA-Z]+)\s*=\s*({.*?}|\'.*?\'|\".*?\"|'
                        '[0-9\s]+\s?|[^=\s]+\s*[0-9]*)')

# Regular expression to extract region type or coordinate system
regex_global = re.compile('^#? *(-?)([a-zA-Z0-9]+)')

# Regular expression to strip parenthesis
regex_paren = re.compile('[()]')

# Regular expression to split coordinate strings
regex_splitter = re.compile('[, ]')


shape_to_region = {}
shape_to_region['pixel'] = {'circle': CirclePixelRegion,
                            'ellipse': EllipsePixelRegion,
                            'rectangle': RectanglePixelRegion,
                            'polygon': PolygonPixelRegion,
                            'circleannulus': CircleAnnulusPixelRegion,
                            'ellipseannulus': EllipseAnnulusPixelRegion,
                            'rectangleannulus': RectangleAnnulusPixelRegion,
                            'line': LinePixelRegion,
                            'point': PointPixelRegion,
                            'text': TextPixelRegion}
shape_to_region['sky'] = {'circle': CircleSkyRegion,
                          'ellipse': EllipseSkyRegion,
                          'rectangle': RectangleSkyRegion,
                          'polygon': PolygonSkyRegion,
                          'circleannulus': CircleAnnulusSkyRegion,
                          'ellipseannulus': EllipseAnnulusSkyRegion,
                          'rectangleannulus': RectangleAnnulusSkyRegion,
                          'line': LineSkyRegion,
                          'point': PointSkyRegion,
                          'text': TextSkyRegion}


@dataclass
class RegionData:
    frame: str
    region_type: str
    shape: str
    shape_params: str
    meta: dict
    visual: dict


@RegionsRegistry.register(Regions, 'read', 'ds9')
def _read_ds9(filename, errors='strict', cache=False):
    """
    Read a DS9 region file in as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.DS9RegionParserError`. 'warn' will raise a
        `~regions.DS9RegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    with get_readable_fileobj(filename, cache=cache) as fh:
        region_string = fh.read()
        return _parse_ds9(region_string, errors=errors)


#def _zz_split_lines(region_str):
#    for line in region_str.split('\n'):
#        # split on all semicolons,
#        # except those between {} braces ({} contain ds9 text strings)
#        for line_ in _split_semicolon(line):
#            yield line_


def _split_lines(region_str):
    lines = []
    for line in region_str.split('\n'):
        # split on all semicolons,
        # except those between {} braces ({} contain ds9 text strings)
        for line_ in _split_semicolon(line):
            lines.append(line_.strip().lower())
    return lines


# NOTES:
#   * wcs[a-z] maps to icrs (print warning)
# Unsupported:
#   * linear
#   * amplifier
#   * detector
#   * physical

# linear,
coordinate_frames = ['image', 'icrs', 'fk5', 'j2000', 'fk4', 'b1950',
                     'galactic', 'ecliptic']
wcs_frames = ['wcs', 'wcs0'] + [f'wcs{letter}'
                                for letter in string.ascii_lowercase]
coordinate_frames += wcs_frames

regex_frame_or_shape = re.compile('^#? *([+-]?)([a-zA-Z0-9]+)')

# DS9 language specification. This defines how a certain region is read.
# language_spec = {'point': (coordinate, coordinate),
#                     'text': (coordinate, coordinate),
#                     'circle': (coordinate, coordinate, radius),
#                     # This is a special case to deal with n elliptical annuli
#                     'ellipse': itertools.chain((coordinate, coordinate),
#                                             itertools.cycle((radius,))),
#                     'box': (coordinate, coordinate, width, height, angle),
#                     'polygon': itertools.cycle((coordinate,)),
#                     'line': (coordinate, coordinate, coordinate, coordinate),
#                     'annulus': itertools.chain((coordinate, coordinate),
#                                             itertools.cycle((radius,)))}


ds9_shapes = ['circle', 'ellipse', 'box', 'annulus', 'polygon', 'line',
              'point', 'text', 'composite']

shape_templates = {'circle': ('circle',
                                '{center},{radius}'),
                    'ellipse': ('ellipse',
                                '{center},{width},{height}'
                                ',{angle}'),
                    'rectangle': ('box',
                                    '{center},{width},{height},{angle}'),
                    'circleannulus': ('annulus',
                                        '{center},{inner_radius},'
                                        '{outer_radius}'),
                    'ellipseannulus': ('ellipse',
                                        '{center},{inner_width},'
                                        '{inner_height},'
                                        '{outer_width},'
                                        '{outer_height},{angle}'),
                    'rectangleannulus': ('box',
                                        '{center},{inner_width},'
                                        '{inner_height},{outer_width},'
                                        '{outer_height},{angle}'),
                    'polygon': ('polygon',
                                '{vertices}'),
                    'line': ('line',
                            '{start},{end}'),
                    'point': ('point', '{center}'),
                    'text': ('text', '{center}')}


@RegionsRegistry.register(Regions, 'parse', 'ds9')
def _parse_ds9(region_str, errors=None):
    region_data = _parse_region_data(region_str)

    #regions = []
    #for region_data_ in region_data:
    #    regions.extend(_make_region(region_data_))
    #return regions

    return region_data


def _parse_region_data(region_str):
    global_meta = {}
    frame = None
    region_data = []
    composite_meta = ''

    allowed_frames_shapes = coordinate_frames + ds9_shapes

    for line in _split_lines(region_str):
        # skip blank lines
        if not line:
            continue

        # skip comments
        if (line.startswith('#')
                and not line.startswith(('# text(', '# composite('))):
            continue

        # ds9 region files can have multiple (including successive)
        # global lines
        if line.startswith('global'):
            global_meta.update(_parse_meta(line[7:]))
            continue

        match = regex_frame_or_shape.search(line)
        if match is None:
            raise ValueError(f'Unable to parse line "{line}".')
        include_symbol = match.groups()[0]
        frame_or_shape = match.groups()[1]

        if frame_or_shape not in allowed_frames_shapes:
            raise ValueError(f'Unable to parse line "{line}".')

        if frame_or_shape in coordinate_frames:
            frame = frame_or_shape
            continue

        if frame_or_shape in ds9_shapes:
            shape = frame_or_shape
            if frame is None:
                raise ValueError('Coordinate frame was not found for region '
                                 f'"{line}".')

            if frame in wcs_frames:
                warnings.warn('DS9 wcs coordinate frames are unsupported, '
                              'skipping region.')
                continue

            if shape == 'composite':
                idx = line.find('||')
                if idx == -1:
                    raise ValueError(f'unable to parse line "{line}"')
                # composite meta applies to all regions within the
                # composite shape; this will always contain at least
                # "composite=1"
                composite_meta = _parse_meta(line[idx + 2:].strip())

            if include_symbol == '-':
                include = 0
            else:  # '+' or ''
                include = 1
            include_meta = {'include': include}

            params_str, meta_str = _parse_shape(shape, match.span(), line)

            meta, visual = _make_metadata(shape, global_meta, composite_meta,
                                          include_meta, meta_str)

            region_type = 'sky'
            if frame == 'image':
                region_type = 'pixel'

            region_data.append(RegionData(frame, region_type, shape,
                                          params_str, meta, visual))

            # reset composite metadata after the composite region ends
            if '||' not in line and composite_meta:
                composite_meta = {}

    return region_data


def _parse_meta(line):
    #meta_regex = re.compile('([a-zA-Z]+)\s*=\s*({.*?}|\'.*?\'|\".*?\"|'
    #                        '[0-9\s]+\s?|[^=\s]+\s*[0-9]*)')
    meta = {}
    for key, val in meta_regex.findall(line):
        val = val.strip().strip("'").strip('"').lstrip('{').rstrip('}')
        if key not in meta:
            if key == 'tag':
                val = [val]
            meta[key] = val
        else:
            if key == 'tag':
                meta[key].append(val)
            else:
                warnings.warn(f'Found duplicate metadata for "{key}", '
                              'skipping')
    return meta


def _make_metadata(shape, global_meta, composite_meta, include_meta, meta_str):
    all_meta = global_meta.copy()
    all_meta.update(composite_meta)
    all_meta.update(include_meta)
    all_meta.update(_parse_meta(meta_str))

    unsupported = ('line', 'ruler')
    # TODO: support text annotations for all DS9 regions
    #ds9_visual_keys = ('color', 'dashlist', 'dash', 'width', 'font', 'fill',
    #                   'point', 'text')
    ds9_visual_keys = ('color', 'dashlist', 'dash', 'width', 'font', 'fill',
                       'point')

    meta = {}
    visual = {}
    for key, value in all_meta.items():
        if key in unsupported:
            if key == 'line' and '1' not in value:  # ignore this special case
                continue
            warnings.warn(f'DS9 meta "{key}={value}" is unsupported and '
                          'will be dropped')

        try:
            value = int(value)
        except (ValueError, TypeError):
            pass

        if key in ds9_visual_keys:
            visual[key] = value
        else:
            meta[key] = value

    visual = _translate_visual_metadata(visual, shape)

    return meta, visual


def _translate_visual_metadata(visual_meta, shape):
    """
    Translate ds9 visual metadata to mpl.
    """
    meta = visual_meta.copy()

    dash = meta.pop('dash', 0)
    dashlist = meta.pop('dashlist', None)
    if int(dash) == 1:
        if shape == 'point':
            warnings.warn('dashed lines are unsupported for DS9 point '
                          'regions')

        meta['linestyle'] = 'dashed'
        if dashlist is not None:
            meta['dashes'] = [int(i) for i in dashlist.split()]

    # "point=symbol [size]"; size is optional, symbol is not
    point = meta.pop('point', None)
    if point is not None:
        point_ = point.split()
        if len(point_) == 1:
            meta['marker'] = point_[0]
        elif len(point_) == 2:
            meta['marker'], meta['markersize'] = point_
        else:
            raise ValueError(f'invalid point data "{point}"')

    font = meta.pop('font', None)
    if font is not None:
        (meta['fontname'], meta['fontsize'], meta['fontweight'],
         meta['fontstyle']) = font.split()
        meta['fontsize'] = int(meta['fontsize'])
        meta['fontstyle'].replace('roman', 'normal')

    if shape == 'point':
        width = meta.pop('width', None)
        if width is not None:
            meta['markeredgewidth'] = width

    return meta


def _parse_shape(shape, span, line):
    full_line = line
    line = full_line[span[1]:]

    # ds9 writes out text regions in this odd format
    if shape == 'text' and full_line.startswith('# text'):
        idx = line.find(' ')
        if idx == -1:
            raise ValueError(f'unable to parse line "{line}"')
        meta_str = line[idx + 1:]
        params_str = line[:idx]

    else:
        hash_idx = line.find('#')
        if hash_idx == -1:
            params_str = line
            meta_str = ''  # no metadata found
        else:
            params_str = line[:hash_idx]
            meta_str = line[hash_idx + 1:]

    params_str = params_str.strip(' |')  # trailing space and | chars
    params_str = re.sub('[()]', '', params_str).lower()
    meta_str = meta_str.strip()

    return params_str, meta_str





template = {'point': ('coord', 'coord'),
            'text': ('coord', 'coord'),
            'circle': ('coord', 'coord', 'length'),
            #'ellipse': ('coord', 'coord', 'length', 'length', 'angle'),
            #'box': ('coord', 'coord', 'length', 'length', 'angle'),
            'ellipse': ('coord', 'coord', itertools.cycle(('length',))),
            'box': ('coord', 'coord', itertools.cycle(('length',))),
            'polygon': itertools.cycle(('coord',)),
            'line': ('coord', 'coord', 'coord', 'coord'),
            'annulus': itertools.chain(('coord', 'coord'),
                                       itertools.cycle(('length',)))}




def _parse_shape_params(region_data):

    region_type = region_data.region_type
    shape = region_data.shape
    shape = 'ellipse'
    params = [val for val in re.split(r'\s|\,', region_data.shape_params)
              if val]

    nparams = len(params)
    nshapes = 1
    if shape in ('ellipse', 'box') and nparams > 5:
        if nparams % 2 != 1:
            raise ValueError(f'incorrect number of parameters ({nparams}) '
                             'for shape "{shape}"')
        nshapes += (nparams - 5) // 2

    shape_params = []
    for idx, (param_type, value) in enumerate(zip(template[shape], params)):
        if nshapes > 1 and idx == nparams - 1:
            param_type = 'angle'
        shape_params.append((param_type, value))

    # create multiple shapes for multi-ellipse or multi-box regions
    if nshapes > 1:
        tmp_params = []
        for i in range(nshapes):
            idx = (i + 1) * 2
            params = (shape_params[0], shape_params[1], shape_params[idx],
                      shape_params[idx + 1], shape_params[-1])
            tmp_params.append(params)
        shape_params = tmp_params

    return shape_params



def _make_region(region_data):
    region_type = region_data.region_type
    shape = region_data.shape
    shape_params_list = _parse_shape_params(region_data)

    regions = []
    for shape_params in shape_params_list:
        region = shape_to_region[region_type][shape](*shape_params)
        region.meta = RegionMeta(region_data.meta)
        region.visual = RegionVisual(region_data.visual)
        regions.append(region)

    return regions


def _find_text_delim_idx(region_str):
    """
    Find the indices of the DS9 text field delimiters ({}, '', or "") in
    a string.
    """
    pattern = re.compile(r'(text\s*=\s*[{\'"])')
    idx0 = []
    delim = []
    start_idx = []
    for match in pattern.finditer(region_str):
        idx0.append(match.start())
        end_delim = match.group()[-1]
        if end_delim == '{':
            end_delim = '}'
        delim.append(end_delim)
        start_idx.append(match.span()[1])

    idx1 = []
    for sidx, char in zip(start_idx, delim):
        idx1.append(region_str.find(char, sidx))

    return idx0, idx1


def _split_semicolon(region_str):
    r"""
    Split a DS9 region string on semicolons.

    The line is not split on semicolons found in a text field.

    This turned out to be a very tricky problem (attempts with regex
    failed):

      * text fields are delimited by {}, '', or ""
      * the text delimiters do not have to be consistent within a file
      * region strings can have unpaired ' (arcmin) or " (arcsec)
      * region strings can have "#" in a color value (cannot split on #)
      * the text field can contain "text={str}", "text='str'",
        'test="str"' etc., e.g., text={text="hello"}
      * the delimiter characters (escaped or not) used by the text field
        are not allowed within the text field, e.g., "text={my
        field\{test\}}" and "text='my field \'test\''" are invalid.
        However, "text={my field, 'test'}" is valid

    This code finds the text field delimiters and then finds the indices
    of the opening and closing delimiters. Text fields that contain
    "text={str}", etc. are included (as a smaller range between the
    larger text field range), but this is fine because all we need are
    index ranges where to exclude semicolons (for splitting). Semicolons
    found at indices between the open/close delimiter indices are
    excluded from splitting.
    """
    idx0, idx1 = _find_text_delim_idx(region_str)

    semi_idx = [pos for pos, char in enumerate(region_str) if char == ';']
    fidx = []
    for i in semi_idx:
        for i0, i1 in zip(idx0, idx1):
            if i0 <= i <= i1:
                break
        else:
            fidx.append(i + 1)
    fidx.insert(0, 0)

    return [region_str[i:j].rstrip(';')
            for i, j in zip(fidx, fidx[1:] + [None])]
