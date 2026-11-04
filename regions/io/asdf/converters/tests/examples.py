# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from regions.core import PixCoord

parameters = {
    'PixCoord': ['x', 'y'],
    'CirclePixelRegion': ['center', 'radius', 'meta', 'visual'],
    'CircleSkyRegion': ['center', 'radius', 'meta', 'visual'],
    'CircleSphericalSkyRegion': ['center', 'radius', 'meta', 'visual'],
    'RegionMeta': ['background', 'comment', 'component', 'composite', 'corr',
                   'delete', 'edit', 'fixed', 'frame', 'highlite', 'include',
                   'label', 'line', 'move', 'name', 'range', 'restfreq',
                   'rotate', 'select', 'source', 'tag', 'text', 'textrotate',
                   'type', 'veltype'],
    'RegionVisual': ['color', 'dash', 'dashlist', 'default_style', 'edgecolor',
                     'facecolor', 'fill', 'font', 'fontname', 'fontsize',
                     'fontstyle', 'fontweight', 'labelcolor', 'labeloff',
                     'labelpos', 'line', 'linestyle', 'linewidth', 'marker',
                     'markeredgewidth', 'markersize', 'rotation', 'symbol',
                     'symsize', 'symthick', 'textangle', 'textrotate',
                     'usetex',
                     ],
    'EllipsePixelRegion': ['center', 'width', 'height', 'meta', 'visual'],
    'EllipseSkyRegion': ['center', 'width', 'height', 'meta', 'visual'],
    'RectanglePixelRegion': ['center', 'width', 'height', 'meta', 'visual'],
    'RectangleSkyRegion': ['center', 'width', 'height', 'meta', 'visual'],
}


def pixcoord():
    return PixCoord(1, 2), parameters['PixCoord']


def pixcoord_array():
    return PixCoord(np.array([1, 2]), np.array([3, 4])), parameters['PixCoord']


def circle_pixel_region():
    from regions.shapes.circle import CirclePixelRegion
    return (CirclePixelRegion(PixCoord(1, 2), 3),
            parameters['CirclePixelRegion'])


def circle_sky_region():
    from regions.shapes.circle import CircleSkyRegion
    center = SkyCoord(ra=1 * u.deg, dec=2 * u.deg, frame='icrs')
    return (CircleSkyRegion(center, 3 * u.arcsec),
            parameters['CircleSkyRegion'])


def circle_spherical_sky_region():
    from regions.shapes.circle import CircleSphericalSkyRegion
    center = SkyCoord(ra=1 * u.deg, dec=2 * u.deg, frame='icrs')
    return (CircleSphericalSkyRegion(center, 3 * u.arcsec),
            parameters['CircleSphericalSkyRegion'])


def meta():
    from regions.core.metadata import RegionMeta
    meta = {
        'background': True,
        'comment': 'This is a test',
        'corr': 'LL',
        'delete': True,
        'edit': True,
        'fixed': True,
        'frame': 'LSRD',
        'include': True,
        'label': 'Test',
        'line': [1, 1],
        'move': True,
        'range': [-320 * u.m / u.s, -330 * u.m / u.s],
        'restfreq': 1.42 * u.GHz,
        'rotate': True,
        'select': True,
        'source': False,
        'tags': '',
        'veltype': 'RADIO',
    }
    return RegionMeta(**meta), parameters['RegionMeta']


def visual_meta():
    from regions.core.metadata import RegionVisual
    visual_meta = {
        'color': 'red',
        'linestyle': '-',
        'linewidth': 1.0,
        'fill': True,
    }
    return RegionVisual(**visual_meta), parameters['RegionVisual']


def ellipse_pixel_region():
    from regions.shapes.ellipse import EllipsePixelRegion
    return (EllipsePixelRegion(PixCoord(1, 2), 3, 4),
            parameters['EllipsePixelRegion'])


def ellipse_sky_region():
    from regions.shapes.ellipse import EllipseSkyRegion
    center = SkyCoord(ra=1 * u.deg, dec=2 * u.deg, frame='icrs')
    return (EllipseSkyRegion(center, 3 * u.arcsec, 4 * u.arcsec),
            parameters['EllipseSkyRegion'])


def rectangle_pixel_region():
    from regions.shapes.rectangle import RectanglePixelRegion
    return (RectanglePixelRegion(PixCoord(1, 2), 3, 4),
            parameters['RectanglePixelRegion'])


def rectangle_sky_region():
    from regions.shapes.rectangle import RectangleSkyRegion
    center = SkyCoord(ra=1 * u.deg, dec=2 * u.deg, frame='icrs')
    return (RectangleSkyRegion(center, 3 * u.arcsec, 4 * u.arcsec),
            parameters['RectangleSkyRegion'])
