# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.exceptions import AstropyUserWarning
import numpy as np

try:
    import matplotlib  # noqa
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

__all__ = ['DS9RegionParserWarning', 'DS9RegionParserError']


class DS9RegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for DS9 region parsing.
    """


class DS9RegionParserError(ValueError):
    """
    A generic error class for DS9 region parsing.
    """


if not HAS_MATPLOTLIB:
    boxcircle = '8'  # octagon
else:
    import matplotlib.path as mpath

    vertices = np.array([[0., -1.], [0.2652031 , -1.],
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
                         [-1., 1.], [-1., -1.], [0., -1.]])
    codes = np.array([1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                      4, 4, 4, 4, 4, 4, 4, 79, 1, 2, 2, 2, 2, 79],
                     dtype=np.uint8)
    boxcircle = mpath.Path(vertices, codes)

# mapping to matplotlib marker symbols, also compatible with CRTF.
valid_symbols_ds9 = {'circle': 'o',
                     'box': 's',
                     'diamond': 'D',
                     'x': 'x',
                     'cross': '+',
                     'arrow': '^',
                     'boxcircle': boxcircle}
