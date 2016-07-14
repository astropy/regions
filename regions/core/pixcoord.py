# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numbers
import numpy as np
from astropy.coordinates import SkyCoord

__all__ = ['PixCoord']


class PixCoord(object):
    """Pixel coordinates.

    This class can represent scalar or array pixel coordinates.

    TODO: see examples below and / or in high-level docs section.

    The data members are either numbers or Numpy arrays
    (not `~astropy.units.Quantity` objects with unit "pixel").

    Given a `astropy.wcs.WCS` object, it can be transformed to and from a
    `~astropy.coordinates.SkyCoord` object.

    TODO: 2-dim or n-dim arrays of pixel coordinates are not supported yet.
    Having those would be very useful, e.g. for arrays of pixel coordinates in
    a 2-dim image.

    Parameters
    ----------
    x : float or array-like
        Pixel coordinate x value
    y : float or array-like
        Pixel coordinate y value
    origin : int
        Origin of pixel coordinates (usually 0 or 1)
    mode : {'all', 'wcs'}
        Whether to do the transformation including distortions ('all')
        or only including only the core WCS transformation ('wcs').

    Examples
    --------
    >>> from regions import PixCoord

    Scalar pixel coordinate:

    >>> pixcoord = PixCoord(x=1, y=11)
    >>> print(pixcoord)
    PixCoord
    x : 1
    y : 11

    Array pixel coordinates:

    >>> pixcoord = PixCoord(x=[1, 2], y=[11, 12])
    >>> print(pixcoord)
    PixCoord
    x : [1 2]
    y : [11 12]
    """
    # We define defaults origin and mode here for WCS transforms.
    # These defaults are used in a few places as default keyword arguments to methods.
    # The purpose is to ensure consistency across the codebase.
    # They should not be modified by users, it's an implementation detail.
    # TODO: is this the best place to put those defaults?
    # Should we put them somewhere else?
    # Or remove them and try to achieve consistency manually?
    _DEFAULT_ORIGIN = 0
    _DEFAULT_MODE = 'all'

    def __init__(self, x, y, origin=_DEFAULT_ORIGIN, mode=_DEFAULT_MODE):
        self.x = self._standardize_coordinate(x)
        self.y = self._standardize_coordinate(y)
        self.origin = origin
        self.mode = mode

    @staticmethod
    def _standardize_coordinate(val):
        """Standardize coordinate.

        Output should be either:
        * Scalar Python number
        * Numpy array

        That's it, nothing else allowed!
        """
        if isinstance(val, numbers.Number):
            return val
        else:
            return np.array(val)

    @property
    def isscalar(self):
        """Is this pixcoord a scalar? (a bool property)"""
        return np.isscalar(self.x)

    def __repr__(self):
        d = dict(name=self.__class__.__name__, x=self.x, y=self.y)
        s = '{name}\nx : {x}\ny : {y}'
        return s.format(**d)

    def __len__(self):
        """Define len(pixcoord) for array-valued pixcoord"""
        return len(self.x)

    def __iter__(self):
        """Allows iteration for array-valued pixcoord

        Yields scalar `PixCoord` objects.
        """
        for (x, y) in zip(self.x, self.y):
            yield PixCoord(x=x, y=y)

    def __getitem__(self, key):
        """Define indexing and slicing."""
        if self.isscalar:
            raise IndexError('Scalar PixCoord cannot be indexed or sliced.')

        # Let Numpy do the slicing 
        x = self.x[key]
        y = self.y[key]
        return PixCoord(x=x, y=y)

    def to_sky(self, wcs):
        """Convert this `PixCoord` to `~astropy.coordinates.SkyCoord`.

        Calls :meth:`astropy.coordinates.SkyCoord.from_pixel`.
        """
        return SkyCoord.from_pixel(
            xp=self.x, yp=self.y, wcs=wcs,
            origin=self.origin, mode=self.mode,
        )

    @classmethod
    def from_sky(cls, skycoord, wcs, origin=_DEFAULT_ORIGIN, mode=_DEFAULT_MODE):
        """Create `PixCoord` from `~astropy.coordinates.SkyCoord`.

        Calls :meth:`astropy.coordinates.SkyCoord.to_pixel`.
        """
        x, y = skycoord.to_pixel(wcs=wcs, origin=origin, mode=mode)
        return cls(x=x, y=y, origin=origin, mode=mode)

    def to_shapely(self):
        """Convert this coord object to a `shapely.geometry.Point` object.

        The ``origin`` and ``mode`` attributes are discarded,
        i.e. this is a lossy conversion.
        """
        from shapely.geometry import Point
        return Point(self.x, self.y)

    @classmethod
    def from_shapely(cls, point, origin=_DEFAULT_ORIGIN, mode=_DEFAULT_MODE):
        return cls(x=point.x, y=point.y, origin=origin, mode=mode)
