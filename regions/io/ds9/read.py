# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import deepcopy
from dataclasses import dataclass
import itertools
import re
import string
import warnings

from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
from astropy.utils.data import get_readable_fileobj
from astropy.utils.exceptions import AstropyUserWarning

from ...core import Regions, RegionMeta, RegionVisual, PixCoord
from ...core.registry import RegionsRegistry
from .core import ds9_valid_symbols

from ...shapes import (CirclePixelRegion, CircleSkyRegion,
                       EllipsePixelRegion, EllipseSkyRegion,
                       RectanglePixelRegion, RectangleSkyRegion,
                       PolygonPixelRegion, PolygonSkyRegion,
                       CircleAnnulusPixelRegion, CircleAnnulusSkyRegion,
                       EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion,
                       RectangleAnnulusPixelRegion, RectangleAnnulusSkyRegion,
                       LinePixelRegion, LineSkyRegion,
                       PointPixelRegion, PointSkyRegion,
                       TextPixelRegion, TextSkyRegion)

__all__ = []


@RegionsRegistry.register(Regions, 'read', 'ds9')
def _read_ds9(filename, cache=False):
    """
    Read a DS9 region file in as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

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
        return _parse_ds9(region_string)


@RegionsRegistry.register(Regions, 'parse', 'ds9')
def _parse_ds9(region_str):
    """
    Parse a DS9 region string.

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    region_data = _parse_region_data(region_str)
    return region_data

    regions = []
    for region_data_ in region_data:
        regions.extend(_make_region(region_data_))
    return regions


@dataclass
class _RegionData:
    """
    Class to hold data used to initialize a Region object.

    Data for multi-annulus regions is stored in a single object.
    """
    frame: str
    region_type: str
    shape: str
    shape_params: str
    meta: dict
    visual: dict


def _split_lines(region_str):
    """
    Split a region string on newlines and all semicolons, except those
    between {} braces ({} contains ds9 text strings).

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    lines : list of str
        A list of strings.
    """
    lines = []
    for line in region_str.split('\n'):
        for line_ in _split_semicolon(line):
            lines.append(line_.strip())
    return lines


def _parse_region_data(region_str):
    """
    Parse region data.

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    lines : list of `_RegionData`
        A list of `_RegionData` objects. Data for multi-annulus regions
        is stored in a single object.
    """
    global_meta = {}
    composite_meta = ''
    frame = None
    region_data = []

    ds9_shapes = ['circle', 'ellipse', 'box', 'annulus', 'polygon', 'line',
                  'point', 'text', 'composite']
    unsupported_shapes = ['vector', 'ruler', 'compass', 'projection',
                          'epanda', 'bpanda']

    coordinate_frames = ['image', 'icrs', 'fk5', 'j2000', 'fk4', 'b1950',
                         'galactic', 'ecliptic']
    unsupported_frames = ['linear', 'amplifier', 'detector', 'physical']
    wcs_frames = ['wcs', 'wcs0'] + [f'wcs{letter}'
                                    for letter in string.ascii_lowercase]
    unsupported_frames += wcs_frames

    allowed_frames_shapes = coordinate_frames + ds9_shapes
    regex_frame_or_shape = re.compile('^#? *([+-]?)([a-zA-Z0-9]+)')

    for line in _split_lines(region_str):  # split on semicolons & newlines
        # skip blank lines
        if not line:
            continue

        # skip comments
        if (line.startswith('#')
                and not line.startswith(('# text(', '# composite('))):
            continue

        original_line = line  # used to parse text and tag fields (keep case)
        line = line.lower()

        # ds9 region files can have multiple (including successive)
        # global lines
        if line.startswith('global'):
            global_meta.update(_parse_metadata(original_line[7:]))
            continue

        match = regex_frame_or_shape.search(line)
        if match is None:
            raise ValueError(f'Unable to parse line "{line}".')
        include_symbol = match.groups()[0]
        frame_or_shape = match.groups()[1]

        if frame_or_shape not in allowed_frames_shapes:
            raise ValueError(f'Unable to parse line "{line}".')

        if frame_or_shape in coordinate_frames:
            # frame persists for subsequent regions until changed
            frame = frame_or_shape
            continue

        if frame_or_shape in ds9_shapes:
            shape = frame_or_shape

            if shape in unsupported_shapes:
                warnings.warn(f'DS9 shape "{shape}" is unsupported, '
                              'skipping the corresponding region.',
                              AstropyUserWarning)
                continue

            # a frame must be defined before any shape(s)
            if frame is None:
                raise ValueError('Coordinate frame was not found for region '
                                 f'"{line}".')

            if frame in unsupported_frames:
                warnings.warn(f'DS9 coordinate frame {frame} is unsupported, '
                              'skipping the corresponding region.',
                              AstropyUserWarning)
                continue

            if shape == 'composite':
                idx = line.find('||')
                if idx == -1:
                    raise ValueError('unable to parse line with composite '
                                     f'shape: "{line}"')
                # composite metadata applies to all regions within the
                # composite shape; this will always contain at least
                # "composite=1"
                composite_meta = _parse_metadata(line[idx + 2:].strip())

            # NOTE: include=1/0 in metadata overrides the leading
            #       "-/+" symbol
            if include_symbol == '-':
                include = 0
            else:  # '+' or ''
                include = 1
            include_meta = {'include': include}

            params_str, meta_str = _parse_shape_line(shape, original_line,
                                                     match.span())

            region_meta = _parse_metadata(meta_str)
            meta, visual = _define_region_metadata(shape, global_meta,
                                                   composite_meta,
                                                   include_meta, region_meta)

            region_type = 'sky'
            if frame == 'image':
                region_type = 'pixel'

            region_data.append(_RegionData(frame, region_type, shape,
                                           params_str, meta, visual))

            # reset composite metadata after the composite region ends
            if '||' not in line and composite_meta:
                composite_meta = {}

    return region_data


def _parse_metadata(metadata_str):
    """
    Parse metadata for a single ds9 region.

    Parameters
    ----------
    metadata_str : str
        The region metadata (e.g., everything after the "#").

    Returns
    -------
    metadata : dict
        The region metadata as a dictionary.
    """
    metadata_regex = re.compile(r'([a-zA-Z]+)\s*=\s*({.*?}|\'.*?\'|\".*?\"|'
                                r'[0-9\s]+\s?|[^=\s]+\s*[0-9]*)')

    metadata = {}
    for key, val in metadata_regex.findall(metadata_str):
        val = val.strip().strip("'").strip('"').lstrip('{').rstrip('}')
        key = key.lower()
        if key not in metadata:
            if key == 'tag':
                val = [val]  # tag value is always a list
            metadata[key] = val
        else:
            if key == 'tag':
                metadata[key].append(val)
            else:
                warnings.warn(f'Found duplicate metadata for "{key}", '
                              'skipping', AstropyUserWarning)
    return metadata


def _define_region_metadata(shape, global_meta, composite_meta, include_meta,
                            region_meta):
    """
    Define the meta and visual dictionaries of metadata for the region.

    Parameters
    ----------
    shape : str
        The DS9 region shape (e.g., 'circle').

    global_meta : dict
        The global metadata.

    composite_meta : dict
        The metadata for a composite region.

    include_meta : dict
        The include/exclude metadata.

    region_meta : dict
        The region metadata.

    Returns
    -------
    meta, visual : tuple of dict
        The meta and visual dictionaries of region metadata. The visual
        metadata is translated to matplotlib keywords.
    """
    all_meta = global_meta.copy()
    all_meta.update(composite_meta)
    all_meta.update(include_meta)
    # region_meta must come after include_meta because include=1/0 in
    # metadata overrides the leading "-/+" include symbol
    all_meta.update(region_meta)

    unsupported = ('line', 'ruler')
    # TODO: include text in visual keys to support text annotations for
    # all DS9 regions?
    ds9_visual_keys = ('color', 'dash', 'dashlist', 'fill', 'font', 'point',
                       'textangle', 'textrotate', 'width')

    meta = {}
    visual = {}
    for key, value in all_meta.items():
        if key in unsupported:
            if key == 'line' and '1' not in value:  # ignore this special case
                continue
            warnings.warn(f'DS9 meta "{key}={value}" is unsupported and '
                          'will be dropped', AstropyUserWarning)

        try:
            value = int(value)
        except (ValueError, TypeError):
            pass

        if key in ds9_visual_keys:
            visual[key] = value
        else:
            meta[key] = value

    visual = _translate_visual_metadata(shape, visual)

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

    dash = meta.pop('dash', 0)
    dashlist = meta.pop('dashlist', None)
    if int(dash) == 1:
        if shape == 'point':
            warnings.warn('dashed lines are unsupported for DS9 point '
                          'regions', AstropyUserWarning)

        if dashlist is not None:
            dashes = tuple(int(i) for i in dashlist.split())
            meta['linestyle'] = (0, dashes)
        else:
            meta['linestyle'] = 'dashed'

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
        meta['fontsize'] = int(meta['fontsize'])
        meta['fontstyle'].replace('roman', 'normal')

    if shape == 'point':
        width = meta.pop('width', None)
        if width is not None:
            meta['markeredgewidth'] = width

    if shape == 'text':
        textangle = meta.pop('textangle', None)
        if textangle is not None:
            meta['rotation'] = textangle

    if shape not in ('point', 'line', 'text'):
        color = meta.pop('color', None)
        if color is not None:
            meta['facecolor'] = color
            meta['edgecolor'] = color

    return meta


def _parse_shape_line(shape, line, span):
    """
    Parse a shape line of a DS9 file.

    Parameters
    ----------
    shape : str
        The DS9 region shape (e.g., 'circle').

    line : str
        A line defining a DS9 region.

    span : 2 tuple
        A tuple containing the (start, end) positions of the match
        defining the shape.

    Returns
    -------
    shape_params_str : str
        The region shape parameters as a string.

    meta_str : str
        The region metadata as a string.
    """
    full_line = line
    line = full_line[span[1]:]  # starts with the shape parameters

    # ds9 writes out text regions in this odd (undocumented) format
    if shape == 'text' and full_line.lower().startswith('# text'):
        idx = line.find(' ')
        if idx == -1:
            raise ValueError(f'unable to parse line "{line}"')
        meta_str = line[idx + 1:]
        shape_params_str = line[:idx]

    else:
        # split line into shape parameters and metadata
        parts = line.split('#', 1)  # color value can contain a #
        if len(parts) == 1:
            shape_params_str = line
            meta_str = ''  # no metadata found
        else:
            shape_params_str, meta_str = parts

    # strip trailing space and | chars
    shape_params_str = shape_params_str.strip(' |')
    shape_params_str = re.sub('[()]', '', shape_params_str).lower()
    meta_str = meta_str.strip()

    return shape_params_str, meta_str


def _parse_pixel_coord(param_str):
    invalid_chars = ('d', 'r', 'p', ':', 'h', 'm', 's')
    for char in invalid_chars:
        if char in param_str:
            warnings.warn('Cannot parse pixel region position coordinates',
                          AstropyUserWarning)
            return None

    if param_str[-1] == 'i':
        param_str = param_str[:-1]

    # DS9 uses 1-indexed pixels
    return float(param_str) - 1


def _parse_sky_coord(param_str, frame, index):
    invalid_chars = ('i', 'p')
    for char in invalid_chars:
        if char in param_str:
            warnings.warn('Cannot parse sky region position coordinates',
                          AstropyUserWarning)
            return None

    if param_str[-1] == 'r':
        return Angle(param_str[:-1], unit=u.radian)

    elif param_str[-1] == 'd':
        return Angle(param_str[:-1], unit=u.degree)

    elif 'd' in param_str or 'h' in param_str:
        return Angle(param_str)

    elif ':' in param_str:
        if index % 2 == 0 and frame not in ('galactic',):
            return Angle(param_str, u.hourangle)
        else:
            return Angle(param_str, u.degree)

    else:  # unit not specified
        return Angle(float(param_str), unit=u.degree)


def _parse_coord(region_type, param_str, frame, index):
    if region_type == 'pixel':
        return _parse_pixel_coord(param_str)
    else:
        return _parse_sky_coord(param_str, frame, index)


def _parse_angle(param_str):
    invalid_chars = ('p', 'i')
    for char in invalid_chars:
        if char in param_str:
            warnings.warn('Cannot parse sky region parameter',
                          AstropyUserWarning)
            return None

    unit_mapping = {'"': u.arcsec,
                    "'": u.arcmin,
                    'd': u.deg,
                    'r': u.rad}
    if param_str[-1] not in string.digits:
        unit = unit_mapping[param_str[-1]]
        return u.Quantity(float(param_str[:-1]), unit=unit)
    else:
        return u.Quantity(float(param_str), unit=u.degree)


def _parse_size(region_type, param_str):
    if region_type == 'pixel':
        invalid_chars = ('"', "'", 'd', 'r', 'p')
        for char in invalid_chars:
            if char in param_str:
                warnings.warn('Cannot parse pixel region parameters',
                              AstropyUserWarning)
                return None

        if param_str[-1] == 'i':
            param_str = param_str[:-1]

        return float(param_str)

    else:
        # size in angular units
        return _parse_angle(param_str)


def _parse_shape_params(region_data):
    """
    Parse the shape parameters for a region line.

    Parameters
    ----------
    region_data : `_RegionData` instance
        A `_RegionData` instance containing the data for a region line.
        Data for multi-annulus regions is stored in a single object.

    Returns
    -------
    shape_params : list of list(s)
        The shape parameters for each region(s). ``shape_params`` is
        a list of lists of shape parameters. A separate shape list is
        returned for multiple annulus/box/ellipse regions. Otherwise,
        the list contains only one list.
    """
    ds9_template = {'point': ('coord', 'coord'),
                    'text': ('coord', 'coord'),
                    'circle': ('coord', 'coord', 'length'),
                    'ellipse': itertools.chain(('coord', 'coord'),
                                               itertools.cycle(('length',))),
                    'box': itertools.chain(('coord', 'coord'),
                                           itertools.cycle(('length',))),
                    'polygon': itertools.cycle(('coord',)),
                    'line': ('coord', 'coord', 'coord', 'coord'),
                    'annulus': itertools.chain(('coord', 'coord'),
                                               itertools.cycle(('length',)))}

    region_type = region_data.region_type
    shape = region_data.shape
    frame = region_data.frame
    params = [val for val in re.split(r'\s|\,', region_data.shape_params)
              if val]  # split values on space or comma

    nparams = len(params)
    n_annulus = 0

    if shape in ('ellipse', 'box') and nparams > 5:
        if nparams % 2 != 1:
            raise ValueError(f'incorrect number of parameters ({nparams}) '
                             f'for shape "{shape}"')
        n_annulus = ((nparams - 3) // 2) - 1

    if shape in ('annulus',):
        n_annulus = nparams - 3

    if shape in ('ellipse', 'box', 'annulus'):
        # deepcopy to "reset" the cycle iterators
        shape_template = deepcopy(ds9_template[shape])
    else:
        shape_template = ds9_template[shape]

    shape_params = []
    for idx, (param_type, value) in enumerate(zip(shape_template, params)):
        if shape in ('ellipse', 'box') and idx == nparams - 1:
            param_type = 'angle'  # last parameter is always an angle

        if param_type == 'coord':
            param = _parse_coord(region_type, value, frame, idx)
        elif param_type in ('length',):
            param = _parse_size(region_type, value)
            if shape == 'ellipse':
                param *= 2.0  # ds9 uses semi-axis lengths
        elif param_type in ('angle',):
            param = _parse_angle(value)
        else:
            raise ValueError('cannot parse shape parameters')

        shape_params.append(param)

    if n_annulus > 1:
        tmp_params = []
        for i in range(n_annulus):
            idx = (i + 1) * 2
            params = [shape_params[0], shape_params[1], shape_params[idx],
                      shape_params[idx + 1], shape_params[idx + 2],
                      shape_params[idx + 3], shape_params[-1]]
            tmp_params.append(params)
        shape_params = tmp_params
    else:
        shape_params = [shape_params]

    if n_annulus > 0:
        if shape == 'ellipse':
            shape = 'ellipse_annulus'
        elif shape == 'box':
            shape = 'rectangle_annulus'
        else:
            raise ValueError('cannot parse shape parameters')

    return shape, shape_params


def _define_pixel_params(shape, shape_params):
    if shape == 'polygon':
        params = [PixCoord(shape_params[0::2], shape_params[1::2])]
    elif shape == 'line':
        params = [PixCoord(*shape_params[0:2]), PixCoord(*shape_params[2:4])]
    elif shape in ('ellipse_annulus', 'rectangle_annulus'):
        size_params = shape_params[2:-1]
        tmp = [size_params[0::2], size_params[1::2]]
        tmp_flat = [item for sublist in tmp for item in sublist]
        params = [PixCoord(*shape_params[0:2]), *tmp_flat, shape_params[-1]]
    else:
        params = ([PixCoord(shape_params[0], shape_params[1])]
                  + shape_params[2:])
    return params


def _define_sky_params(shape, shape_params, frame):
    frame_mapping = {'image': 'image',
                     'icrs': 'icrs',
                     'fk5': 'fk5',
                     'fk4': 'fk4',
                     'galactic': 'galactic',
                     'geocentrictrueecliptic': 'ecliptic'}

    frame_mapping_rev = {val: key for key, val in frame_mapping.items()}
    frame = frame_mapping_rev[frame]

    if shape == 'polygon':
        params = [SkyCoord(shape_params[0::2], shape_params[1::2],
                           frame=frame)]
    elif shape == 'line':
        params = [SkyCoord(*shape_params[0:2], frame=frame),
                  SkyCoord(*shape_params[2:4], frame=frame)]
    else:
        params = ([SkyCoord(shape_params[0], shape_params[1],
                            frame=frame)] + shape_params[2:])
    return params


def _make_region(region_data):
    """
    region_data : `_RegionData` instance
    """
    pixel_map = {'circle': CirclePixelRegion,
                 'ellipse': EllipsePixelRegion,
                 'box': RectanglePixelRegion,
                 'polygon': PolygonPixelRegion,
                 'annulus': CircleAnnulusPixelRegion,
                 'ellipse_annulus': EllipseAnnulusPixelRegion,
                 'rectangle_annulus': RectangleAnnulusPixelRegion,
                 'line': LinePixelRegion,
                 'point': PointPixelRegion,
                 'text': TextPixelRegion}

    sky_map = {'circle': CircleSkyRegion,
               'ellipse': EllipseSkyRegion,
               'box': RectangleSkyRegion,
               'polygon': PolygonSkyRegion,
               'annulus': CircleAnnulusSkyRegion,
               'ellipse_annulus': EllipseAnnulusSkyRegion,
               'rectangle_annulus': RectangleAnnulusSkyRegion,
               'line': LineSkyRegion,
               'point': PointSkyRegion,
               'text': TextSkyRegion}

    shape_to_region = {}
    shape_to_region['pixel'] = pixel_map
    shape_to_region['sky'] = sky_map

    region_type = region_data.region_type
    shape = region_data.shape
    frame = region_data.frame
    shape, shape_params_list = _parse_shape_params(region_data)

    # define the parameters to initalize a Region
    region_params = []
    for shape_params in shape_params_list:
        if region_type == 'pixel':
            region_params.extend([_define_pixel_params(shape, shape_params)])
        else:
            region_params.extend([_define_sky_params(shape, shape_params,
                                                     frame)])

    regions = []
    for shape_params in region_params:
        # for Text region, we need to add meta['text'] to params
        if shape == 'text':
            #print(region_data.meta)
            #print(region_data.meta['text'])
            shape_params.append(region_data.meta['text'])

        #print(region_type, shape, shape_params)
        #print(shape_to_region[region_type][shape])
        print(region_type, shape, shape_params)

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
