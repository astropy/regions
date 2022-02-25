# Licensed under a 3-clause BSD style license - see LICENSE.rst
import itertools

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

from ..._utils.optional_deps import HAS_MATPLOTLIB

__all__ = []


# mappings from DS9 frames to astropy coordinates frames
ds9_frame_map = {'image': 'image',
                 'icrs': 'icrs',
                 'fk5': 'fk5',
                 'j2000': 'fk5',
                 'fk4': 'fk4',
                 'b1950': 'fk4',
                 'galactic': 'galactic',
                 'ecliptic': 'barycentricmeanecliptic'}

# mappings from DS9 shape to region class
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

ds9_shape_to_region = {}
ds9_shape_to_region['pixel'] = pixel_map
ds9_shape_to_region['sky'] = sky_map


ds9_params_template = {'point': ('coord', 'coord'),
                       'text': ('coord', 'coord'),
                       'circle': ('coord', 'coord', 'length'),
                       'ellipse': itertools.chain(
                           ('coord', 'coord'), itertools.cycle(('length',))),
                       'box': itertools.chain(
                           ('coord', 'coord'), itertools.cycle(('length',))),
                       'polygon': itertools.cycle(('coord',)),
                       'line': ('coord', 'coord', 'coord', 'coord'),
                       'annulus': itertools.chain(
                           ('coord', 'coord'), itertools.cycle(('length',)))}


# mapping from regions shapes to ds9 shape formats
# unsupported ds9 shapes:
# vector, ruler, compass, projection, panda, epanda, bpanda, composite
ds9_shape_templates = {'circle': ('circle',
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


if not HAS_MATPLOTLIB:
    boxcircle = '8'  # octagon
    arrow = '^'  # triangle
else:
    import matplotlib.path as mpath

    vertices = [[0., -1.], [0.2652031, -1.],
                [0.51957987, -0.89463369], [0.70710678, -0.70710678],
                [0.89463369, -0.51957987], [1., -0.2652031],
                [1., 0.], [1., 0.2652031], [0.89463369, 0.51957987],
                [0.70710678, 0.70710678], [0.51957987, 0.89463369],
                [0.2652031, 1.], [0., 1.], [-0.2652031, 1.],
                [-0.51957987, 0.89463369], [-0.70710678, 0.70710678],
                [-0.89463369, 0.51957987], [-1., 0.2652031],
                [-1., 0.], [-1., -0.2652031],
                [-0.89463369, -0.51957987],
                [-0.70710678, -0.70710678],
                [-0.51957987, -0.89463369], [-0.2652031, -1.],
                [0., -1.], [0., -1.], [0., -1.], [1., -1.], [1., 1.],
                [-1., 1.], [-1., -1.], [0., -1.]]
    codes = [1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
             4, 4, 4, 4, 4, 4, 4, 79, 1, 2, 2, 2, 2, 79]
    boxcircle = mpath.Path(vertices, codes)

    arrow_verts = [[-1, 0], [0, 0], [-1, 1], [0, 0], [0, 1]]
    arrow = mpath.Path(arrow_verts, codes=None)


# mapping from ds9 point symbols to matplotlib marker symbols
ds9_valid_symbols = {'circle': 'o',
                     'box': 's',
                     'diamond': 'D',
                     'x': 'x',
                     'cross': '+',
                     'arrow': arrow,
                     'boxcircle': boxcircle}


class DS9ParserError(Exception):
    """
    A custom exception for DS9 parsing errors.
    """
