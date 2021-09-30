# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings

from astropy.coordinates import Angle
from astropy.units import Quantity
from astropy.utils.exceptions import AstropyUserWarning

from ...core import Region, Regions, PixelRegion
from ...core.registry import RegionsRegistry
from ..core import _to_shape_list

__all__ = []


@RegionsRegistry.register(Region, 'serialize', 'ds9')
@RegionsRegistry.register(Regions, 'serialize', 'ds9')
def _serialize_ds9(regions, fmt='.6f', radunit='deg'):
    shapelist = _to_shape_list(regions)
    return shapelist.to_ds9(fmt, radunit)


@RegionsRegistry.register(Region, 'write', 'ds9')
@RegionsRegistry.register(Regions, 'write', 'ds9')
def _write_ds9(regions, filename, fmt='.6f', radunit='deg',
               overwrite=False):
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

    output = _serialize_ds9(regions, fmt=fmt, radunit=radunit)
    with open(filename, 'w') as fh:
        fh.write(output)


def _get_frame_name(region, mapping):
    if isinstance(region, PixelRegion):
        frame = 'image'
    else:
        frame = region.center.frame.name

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


def _get_region_center(region, precision=8):
    if isinstance(region, PixelRegion):
        # pixels (TODO: apply precision?)
        center = f'{region.center.x},{region.center.y}'
    else:
        # decimal degrees
        center = region.center.to_string(precision=precision).replace(' ', ',')
    return center


def _get_shape_params(region, template, precision=8):
    param = {}
    for param_name in region._params:
        if param_name == 'center':
            continue
        value = getattr(region, param_name)
        if isinstance(value, Angle):
            value = value.to_string(unit='deg', decimal=True,
                                    precision=precision)
        elif isinstance(value, Quantity):
            value = value.to_string(unit='deg', precision=precision)[:-4]
        else:
            value = f'{value:.{precision}f}'
        param[param_name] = value

    try:
        param_str = template.format(**param)
    except KeyError as err:
        raise ValueError(
            f'unable to get shape parameters for {region!r}') from err

    return param_str


def _get_region_meta(region, valid_keys):
    region_meta = {**region.meta, **region.visual}
    meta = {}
    for key in region_meta:
        if key in valid_keys:
            meta[key] = region_meta[key]
    return meta


def _translate_ds9_meta(meta):
    """
    Translate metadata from other regions or matplotlib to valid ds9
    meta keys.
    """
    if 'text' in meta:
        meta.pop('label', None)
    else:
        if 'label' in meta:
            meta['text'] = meta.pop('label')

    if 'width' in meta:
        meta.pop('linewidth', None)
    else:
        if 'linewidth' in meta:
            meta['width'] = meta.pop('linewidth')

    if 'point' in meta:
        meta.pop('symbol', None)
        meta.pop('symsize', None)
    else:
        # symsize without symbol is ignored
        if 'symbol' in meta:
            point = f'{meta.pop("symbol")}'
            if 'symsize' in meta:
                point += f' {meta.pop("symsize")}'
            meta['point'] = point

    if 'font' in meta:
        meta.pop('fontname', None)
        meta.pop('fontsize', None)
        meta.pop('fontweight', None)
        meta.pop('fontstyle', None)
    else:
        # fontsize, fontweight, and fontstyle without fontname are ignored
        if 'fontname' in meta:
            font = f'{meta.pop("fontname")}'
        if 'fontsize' in meta:
            font += f'{meta.pop("fontsize")}'
        if 'fontweight' in meta:
            font += f'{meta.pop("fontweight")}'
        if 'fontstyle' in meta:
            font += f'{meta.pop("fontstyle")}'

    if 'dash' in meta:
        if (meta['dash'] == 0) or (meta['dash'] == 1 and 'dashlist' in meta):
            meta.pop('linestyle', None)
            meta.pop('dashes', None)
    else:
        # if dash not in meta, dashlist is ignored
        if meta.get('linestyle', 'solid') in ('dashed', '--'):
            meta['dash'] = 1
            if 'dashes' in meta:
                dashes = meta.pop('dashes')
                meta['dashlist'] = f'{dashes[0]} {dashes[1]}'

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
                                  '{radius}'),
                       'ellipse': ('ellipse',
                                   '{width},{height},{angle}'),
                       'rectangle': ('box',
                                     '{width},{height}{angle}'),
                       'circleannulus': ('annulus',
                                         '{inner_radius},{outer_radius}'),
                       'ellipseannulus': ('ellipse',
                                          '{inner_width},{inner_height},'
                                          '{outer_width},{outer_height},'
                                          '{angle}'),
                       'rectangleannulus': ('box',
                                            '{inner_width},{inner_height},'
                                            '{outer_width},{outer_height},'
                                            '{angle}'),
                       'polygon': ('polygon',
                                   '{vertices}'),
                       'line': ('line',
                                '{start},{end}'),
                       'point': ('point', ''),
                       'text': ('text', '')}

    frame = _get_frame_name(region, mapping=frame_mapping)
    shape, region_type = _get_region_shape(region, shape_templates)
    region_center = _get_region_center(region, precision=precision)
    template = shape_templates[shape][1]
    shape_params = _get_shape_params(region, template, precision=precision)
    shape_str = f'{region_type}({region_center},{shape_params})'

    # ds9 meta keys
    meta_keys = ['text', 'select', 'highlite', 'fixed', 'edit', 'move',
                 'rotate', 'delete', 'include', 'tag']
    visual_keys = ['color', 'textangle', 'textrotate', 'dash', 'dashlist',
                   'width', 'font', 'fill', 'point']

    # meta keys that can be mapped
    # TODO: make specific to mpl artist type?
    other_keys = ['label', 'linewidth', 'fontname', 'fontsize', 'fontweight',
                  'fontstyle', 'symbol', 'symsize', 'linewidth', 'linestyle',
                  'dashes']

    valid_keys = meta_keys + visual_keys + other_keys

    valid_meta = _get_region_meta(region, valid_keys)
    meta = _translate_ds9_meta(valid_meta)

    return {'frame': frame, 'shape': shape_str, 'meta': meta}
