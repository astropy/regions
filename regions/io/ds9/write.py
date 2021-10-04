# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings

from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.utils.exceptions import AstropyUserWarning

from ...core import Region, Regions, PixelRegion, PixCoord
from ...core.registry import RegionsRegistry
from ..core import _to_shape_list
from .core import valid_symbols_ds9

__all__ = []


#@RegionsRegistry.register(Region, 'serialize', 'ds9')
#@RegionsRegistry.register(Regions, 'serialize', 'ds9')
#def _serialize_ds9(regions, fmt='.6f', radunit='deg'):
#    shapelist = _to_shape_list(regions)
#    return shapelist.to_ds9(fmt, radunit)


@RegionsRegistry.register(Region, 'write', 'ds9')
@RegionsRegistry.register(Regions, 'write', 'ds9')
def _write_ds9(regions, filename, *, precision=8, overwrite=False):
    """
    Convert a list of `~regions.Region` to a DS9 string and write to a
    file.

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    filename : str
        The filename in which the string is to be written.

    fmt : str, optional
        A python string format defining the output precision. Default is
        '.6f', which is accurate to 0.0036 arcseconds.

    radunit : str, optional
        The unit of the radius.

    overwrite : bool, optional
        If True, overwrite the output file if it exists. Raises an
        `OSError` if False and the output file exists. Default is False.
    """
    if os.path.lexists(filename) and not overwrite:
        raise OSError(f'{filename} already exists')

    output = _serialize_ds9(regions, precision=precision)
    with open(filename, 'w') as fh:
        fh.write(output)


def _get_frame_name(region, mapping):
    if isinstance(region, PixelRegion):
        frame = 'image'
    else:
        if 'center' in region._params:
            frame = region.center.frame.name
        elif 'vertices' in region._params:
            frame = region.vertices.frame.name
        elif 'start' in region._params:
            frame = region.start.frame.name
        else:
            raise ValueError(f'Unable to get shape parameters for {region!r}')

    if frame not in mapping.keys():
        warnings.warn(f'Cannot serialize region with frame={frame}, skipping',
                      AstropyUserWarning)

    return mapping[frame]


def _get_region_shape(region, mapping):
    shape = region.__class__.__name__.lower().replace('skyregion', '')
    shape = shape.replace('pixelregion', '')
    shape = shape.replace('regularpolygon', 'polygon')

    if shape not in mapping.keys():
        warnings.warn(f'Cannot serialize region shape "{shape}", '
                      'skipping', AstropyUserWarning)

    return shape, mapping[shape][0]


def _get_region_params(region, template, precision=8):
    ellipse_axes = ('width', 'height', 'inner_width', 'inner_height',
                    'outer_width', 'outer_height')
    ellipse_names = ('ellipse', 'ellipseannulus')

    param = {}
    for param_name in region._params:
        if param_name in ('text'):
            continue

        if param_name in ellipse_axes and template[0] in ellipse_names:
            is_ellipse = True
        else:
            is_ellipse = False

        value = getattr(region, param_name)

        if isinstance(value, PixCoord):
            # pixels; DS9's origin is (1, 1)
            # TODO: apply precision?
            if value.isscalar:
                value = f'{value.x + 1},{value.y + 1}'
            else:
                value_str = ''
                for val in value:
                    value_str += f'{val.x + 1},{val.y + 1},'
                value = value_str[:-1]

        elif isinstance(value, SkyCoord):
            # polygon region
            if not value.isscalar:
                value = ' '.join(value.to_string(
                    precision=precision)).replace(' ', ',')
            else:
                value = value.to_string(
                    precision=precision).replace(' ', ',')

        elif isinstance(value, Angle):
            if is_ellipse:
                value /= 2.0
            value = value.to_string(unit='deg', decimal=True,
                                    precision=precision)

        elif isinstance(value, Quantity):
            if is_ellipse:
                value /= 2.0
            value = value.to_string(unit='deg', precision=precision)[:-4]

        else:
            value = f'{value:.{precision}f}'

        param[param_name] = value

    try:
        param_str = template[1].format(**param)
    except KeyError as err:
        raise ValueError(
            f'Unable to get shape parameters for {region!r}') from err

    return param_str


def _remove_invalid_keys(region_meta, valid_keys):
    # TODO: instead of new dict, del region_meta in-place?
    meta = {}
    for key in region_meta:
        if key in valid_keys:
            meta[key] = region_meta[key]
    return meta


# TODO: remove me after the parsers are refactored
def _translate_to_mpl_meta(region):
    """
    Translate region visual metadata to valid mpl keys.
    """
    meta = region.visual.copy()

    dash = meta.pop('dash', 0)
    dashlist = meta.pop('dashlist', None)
    if int(dash) == 1:
        meta['linestyle'] = 'dashed'
        if dashlist is not None:
            meta['dashes'] = [int(i) for i in dashlist.split()]

    font = meta.pop('font', None)
    if font is not None:
        meta['fontname'] = font

    symbol = meta.pop('symbol', None)
    if symbol is not None:
        meta['marker'] = symbol
    symsize = meta.pop('symsize', None)
    if symsize is not None:
        meta['markersize'] = symsize

    # TODO: if region is point
    #width = meta.pop('width', None)
    #if width is not None:
    #    meta['markeredgewidth'] = width

    return meta


def _translate_ds9_meta(meta):
    """
    Translate metadata from other regions or matplotlib to valid ds9
    meta keys.
    """
    if 'include' in meta:
        if meta['include'] is True:
            meta['include'] = 1
        else:
            meta['include'] = 0
    else:
        meta['include'] = 1

    if 'text' in meta:
        meta['text'] = f'{{{meta["text"]}}}'

    linewidth = meta.pop('linewidth', None)
    if linewidth is not None:
        meta['width'] = linewidth

    # point
    markeredgewidth = meta.pop('markeredgewidth', None)
    if markeredgewidth is not None:
        meta['width'] = markeredgewidth

    fontname = meta.pop('fontname', None)
    if fontname is not None:
        fontsize = meta.pop('fontsize', 10)  # default 10
        fontweight = meta.pop('fontweight', 'normal')
        fontstyle = meta.pop('fontstyle', 'roman').replace('normal',
                                                           'roman')
        meta['font'] = f'"{fontname} {fontsize} {fontweight} {fontstyle}"'

    linestyle = meta.pop('linestyle', None)
    if linestyle in ('dashed', '--'):
        meta['dash'] = 1
        dashes = meta.pop('dashes', None)
        if dashes is not None:
            meta['dashlist'] = f'{dashes[0]} {dashes[1]}'

    marker = meta.pop('marker', None)
    if marker is not None:
        symbol_map = {y: x for x, y in valid_symbols_ds9.items()}
        markersize = meta.pop('markersize', 11)
        meta['point'] = f'{symbol_map[marker]} {markersize}'

    return meta


def _serialize_region_ds9(region, precision=8):
    # mapping from astropy frames to ds9 frames
    frame_mapping = {'image': 'image',
                     'icrs': 'icrs',
                     'fk5': 'fk5',
                     'fk4': 'fk4',
                     'galactic': 'galactic',
                     'geocentrictrueecliptic': 'ecliptic'}

    # mapping from regions shapes to ds9 shapes
    # unsupported ds9 shapes:
    # vector, ruler, compass, projection, panda, epanda, bpanda, composite
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

    frame = _get_frame_name(region, mapping=frame_mapping)
    shape, region_type = _get_region_shape(region, shape_templates)
    template = shape_templates[shape]
    region_params = _get_region_params(region, template, precision=precision)
    region_str = f'{region_type}({region_params})'

    # ds9 meta keys
    meta_keys = ['text', 'select', 'highlite', 'fixed', 'edit', 'move',
                 'rotate', 'delete', 'include', 'tag']
    visual_keys = ['color', 'textangle', 'textrotate', 'dash', 'dashlist',
                   'width', 'font', 'fill', 'point']
    valid_keys = meta_keys + visual_keys

    # TODO: remove after parsers are refactored
    visual_meta = _translate_to_mpl_meta(region)
    region_meta = {**region.meta, **visual_meta}

    meta = _translate_ds9_meta(region_meta)
    meta = _remove_invalid_keys(meta, valid_keys)

    return {'frame': frame, 'region': region_str, 'meta': meta}


def _make_meta_str(meta):
    metalist = []
    for key, val in meta.items():
        if key == 'tag':
            metalist.append(' '.join([f'tag={val}' for val in meta[key]]))
        else:
            metalist.append(f'{key}={val}')
    return ' '.join(metalist)


@RegionsRegistry.register(Region, 'serialize', 'ds9')
@RegionsRegistry.register(Regions, 'serialize', 'ds9')
def _serialize_ds9(regions, precision=8):
    region_data = []
    for region in regions:
        region_data.append(_serialize_region_ds9(region, precision=precision))

    # extract common meta to global
    # tag cannot be in global
    metalist_notag = []
    for rd in region_data:
        tmp = rd['meta']
        tmp.pop('tag', None)
        metalist_notag.append(tmp)

    output = '# Region file format: DS9 astropy/regions\n'

    global_meta = dict(set.intersection(*[set(d.items())
                                          for d in metalist_notag]))
    if global_meta:
        global_str = f'global {_make_meta_str(global_meta)}'
        output += f'{global_str}\n'

    #metalist = [rd['meta'] for rd in region_data]
    metalist = []
    for rd in region_data:
        tmp = rd['meta']
        for key in global_meta.keys():
            tmp.pop(key)
        metalist.append(tmp)

    # TODO
    # extract common coord frame for consecutive regions

    for data, meta in zip(region_data, metalist):
        meta_str = _make_meta_str(meta)
        output += f'{data["frame"]}; {data["region"]} # {meta_str}\n'

    print(output)
    return output
