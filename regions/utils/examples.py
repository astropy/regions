# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from astropy.wcs import WCS

__all__ = ['make_example_wcs']


def make_example_wcs(projection='AIT'):
    """Make a WCS object for documentation examples and tests.

    Parameters
    ----------
    projection : {'AIT', 'CAR'}
        WCS projection type

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        World coordinate system transformation object.

    Examples
    --------
    Make an example WCS object:

    >>> from regions import make_example_wcs
    >>> wcs = make_example_wcs(projection='AIT')

    WCS parameters::

        CTYPE : 'GLON-AIT'  'GLAT-AIT'
        CRVAL : 0.0  0.0
        CRPIX : 18.0  9.0
        CDELT : 10.0  10.0
    """
    if not projection in {'AIT', 'CAR'}:
        raise ValueError('Invalid projection: {}'.format(projection))

    wcs = WCS(naxis=2)
    wcs.wcs.crval = 0, 0
    wcs.wcs.crpix = 18, 9
    wcs.wcs.cdelt = 10, 10
    wcs.wcs.ctype = 'GLON-' + projection, 'GLAT-' + projection

    # shape = (36, 18) would give an image that covers the whole sky.

    return wcs
