# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings

from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.utils.exceptions import AstropyUserWarning

from ...core import Region, Regions, PixelRegion, PixCoord
from ...core.registry import RegionsRegistry
from .core import ds9_valid_symbols, ds9_shape_templates

__all__ = []


@RegionsRegistry.register(Region, 'serialize', 'ds9')
@RegionsRegistry.register(Regions, 'serialize', 'ds9')
def _serialize_ds9(regions, precision=8):
    region_data = []
    for region in regions:
        region_data.append(_serialize_region_ds9(region, precision=precision))

    # ds9 file header
    output = '# Region file format: DS9 astropy/regions\n'

    # extract common region metadata and place in the global metadata
    all_meta = []
    for region in region_data:
        region_meta = region['meta']
        region_meta.pop('tag', None)  # "tag" cannot be in global metadata
        all_meta.append(region_meta)

    global_meta = dict(set.intersection(*[set(meta_dict.items())
                                          for meta_dict in all_meta]))
    if global_meta:
        output += f'global {_make_meta_str(global_meta)}\n'

    # define region frame and metadata (removing items that are
    # defined in global metadata)
    frames = []
    metadata = []
    for region in region_data:
        frames.append(region['frame'])
        meta_tmp = region['meta']
        for key in global_meta:
            meta_tmp.pop(key)
        metadata.append(meta_tmp)

    # extract the coordinate frame if it is identical for all regions
    # TODO: extract common coord frame(s) for block(s) of consecutive regions
    frames = set(frames)
    global_frame = None
    if len(frames) == 1:
        global_frame = frames.pop()
        output += f'{global_frame}\n'

    # add line for each region
    for region, region_meta in zip(region_data, metadata):
        meta_str = _make_meta_str(region_meta)
        if global_frame is None:
            output += f'{region["frame"]}; '
        output += f'{region["region"]} # {meta_str}\n'

    return output


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
        The output region filename.

    precision : int, optional
        The level of decimal precision given as the number of decimal
        places.

    overwrite : bool, optional
        If `True`, overwrite the output file if it exists. If `False`
        (default) and the output file exists, an `OSError` is raised.
    """
    if os.path.lexists(filename) and not overwrite:
        raise OSError(f'{filename} already exists')

    output = _serialize_ds9(regions, precision=precision)
    with open(filename, 'w') as fh:
        fh.write(output)


def _get_region_shape(region):
    shape = region.__class__.__name__.lower().replace('skyregion', '')
    shape = shape.replace('pixelregion', '')
    shape = shape.replace('regularpolygon', 'polygon')
    return shape


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
            raise ValueError(f'Unable to get coordinate frame for {region!r}')

    if frame not in mapping.keys():
        warnings.warn(f'Cannot serialize region with frame={frame}, skipping',
                      AstropyUserWarning)

    return mapping[frame]


def _remove_invalid_keys(region_meta, valid_keys):
    # TODO: instead of new dict, del region_meta in-place?
    meta = {}
    for key in region_meta:
        if key in valid_keys:
            meta[key] = region_meta[key]
    return meta


def _make_meta_str(meta):
    metalist = []
    for key, val in meta.items():
        if key == 'tag':  # can have multiple tags; value is always a list
            metalist.append(' '.join([f'tag={val}' for val in meta[key]]))
        else:
            metalist.append(f'{key}={val}')
    return ' '.join(metalist)


def _get_region_params(region, shape_template, precision=8):
    ellipse_axes = ('width', 'height', 'inner_width', 'inner_height',
                    'outer_width', 'outer_height')
    ellipse_names = ('ellipse', 'ellipseannulus')

    param = {}
    for param_name in region._params:
        if param_name in ('text'):
            continue

        value = getattr(region, param_name)

        # DS9 ellipse parameters are serialized as semi-axis lengths, but
        # ellipse region is defined by full axis lengths
        is_ellipse = (shape_template[0] in ellipse_names
                      and param_name in ellipse_axes)
        if not isinstance(value, (PixCoord, SkyCoord)) and is_ellipse:
            # deepcopy to prevent changing value in memory
            value = deepcopy(value) / 2.0  # semi-axis lengths

        if isinstance(value, PixCoord):
            # pixels; ds9's origin is (1, 1)
            # TODO: apply precision to pixel coordinates?
            if value.isscalar:
                value = f'{value.x + 1},{value.y + 1}'
            else:
                value_str = ''
                for val in value:
                    value_str += f'{val.x + 1},{val.y + 1},'
                value = value_str[:-1]

        elif isinstance(value, SkyCoord):
            value = value.to_string(precision=precision)
            # polygon region has multiple SkyCoord
            if not value.isscalar:
                value = ' '.join(value)
            value = value.replace(' ', ',')

        elif isinstance(value, Angle):
            value = value.to_string(unit='deg', decimal=True,
                                    precision=precision)

        elif isinstance(value, Quantity):
            # [:-4] to trim ' deg' from string end
            value = value.to_string(unit='deg', precision=precision)[:-4]

        else:
            value = f'{value:.{precision}f}'

        param[param_name] = value

    try:
        param_str = shape_template[1].format(**param)
    except KeyError as err:
        raise ValueError(
            f'Unable to get shape parameters for {region!r}') from err

    return param_str


def _translate_ds9_meta(region, shape):
    """
    Translate region metadata to valid ds9 meta keys.
    """
    meta = {**region.meta, **region.visual}

    if 'annulus' in shape:
        # ds9 does not allow fill for annulus regions
        meta.pop('fill', None)

    if 'text' in meta:
        meta['text'] = f'{{{meta["text"]}}}'

    linewidth = meta.pop('linewidth', None)
    if linewidth is not None:
        meta['width'] = linewidth

    # point region marker width
    markeredgewidth = meta.pop('markeredgewidth', None)
    if markeredgewidth is not None:
        meta['width'] = markeredgewidth

    marker = meta.pop('marker', None)
    if marker is not None:
        symbol_map = {y: x for x, y in ds9_valid_symbols.items()}
        markersize = meta.pop('markersize', 11)  # default 11
        meta['point'] = f'{symbol_map[marker]} {markersize}'

    fontname = meta.pop('fontname', None)
    if fontname is not None:
        fontsize = meta.pop('fontsize', 10)  # default 10
        fontweight = meta.pop('fontweight', 'normal')  # default normal
        # default roman
        fontstyle = meta.pop('fontstyle', 'roman').replace('normal', 'roman')
        meta['font'] = f'"{fontname} {fontsize} {fontweight} {fontstyle}"'

    linestyle = meta.pop('linestyle', None)
    if linestyle is not None:
        meta['dash'] = 1
    # if linestyle in ('dashed', '--'):
    if isinstance(linestyle, tuple):
        meta['dashlist'] = f'{linestyle[1][0]} {linestyle[1][1]}'

        # dashes = meta.pop('dashes', None)
        # if dashes is not None:
        #     meta['dashlist'] = f'{dashes[0]} {dashes[1]}'

    rotation = meta.pop('rotation', None)
    if rotation is not None:
        meta['textangle'] = rotation
        # meta['textrotate'] = 1

    # ds9 meta keys
    meta_keys = ['background', 'delete', 'edit', 'fixed', 'highlite',
                 'include', 'move', 'rotate', 'select', 'source', 'tag',
                 'text']
    visual_keys = ['color', 'dash', 'dashlist', 'fill', 'font', 'point',
                   'textangle', 'textrotate', 'width']
    valid_keys = meta_keys + visual_keys
    meta = _remove_invalid_keys(meta, valid_keys)

    return meta


def _serialize_region_ds9(region, precision=8):
    # mapping from astropy frames to ds9 frames
    frame_mapping = {'image': 'image',
                     'icrs': 'icrs',
                     'fk5': 'fk5',
                     'fk4': 'fk4',
                     'galactic': 'galactic',
                     'geocentrictrueecliptic': 'ecliptic'}

    frame = _get_frame_name(region, mapping=frame_mapping)

    shape = _get_region_shape(region)
    if shape not in ds9_shape_templates:
        warnings.warn(f'Cannot serialize region shape "{shape}", '
                      'skipping', AstropyUserWarning)

    region_params = _get_region_params(region,
                                       ds9_shape_templates[shape],
                                       precision=precision)

    region_type = ds9_shape_templates[shape][0]
    region_str = f'{region_type}({region_params})'

    region_meta = _translate_ds9_meta(region, shape)

    return {'frame': frame, 'region': region_str, 'meta': region_meta}
