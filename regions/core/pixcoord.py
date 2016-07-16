# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numbers
import numpy as np
from astropy.coordinates import SkyCoord
from .core import _DEFAULT_WCS_ORIGIN, _DEFAULT_WCS_MODE

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

    def __init__(self, x, y):
        self.x = self._standardize_coordinate(x)
        self.y = self._standardize_coordinate(y)

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

    def to_sky(self, wcs, origin=_DEFAULT_WCS_ORIGIN, mode=_DEFAULT_WCS_MODE):
        """Convert this `PixCoord` to `~astropy.coordinates.SkyCoord`.

        Calls :meth:`astropy.coordinates.SkyCoord.from_pixel`.
        See parameter description there.
        """
        return SkyCoord.from_pixel(
            xp=self.x, yp=self.y, wcs=wcs,
            origin=origin, mode=mode,
        )

    @classmethod
    def from_sky(cls, skycoord, wcs, origin=_DEFAULT_WCS_ORIGIN, mode=_DEFAULT_WCS_MODE):
        """Create `PixCoord` from `~astropy.coordinates.SkyCoord`.

        Calls :meth:`astropy.coordinates.SkyCoord.to_pixel`.
        See parameter description there.
        """
        x, y = skycoord.to_pixel(wcs=wcs, origin=origin, mode=mode)
        return cls(x=x, y=y)

    def to_shapely(self):
        """Convert this coord object to a `shapely.geometry.Point` object.
        """
        if not self.isscalar:
            raise TypeError('Scalar PixCoord cannot be converted to Shapely.')

        from shapely.geometry import Point
        return Point(self.x, self.y)

    @classmethod
    def from_shapely(cls, point):
        """Create `PixCoord` from `shapely.geometry.Point` object.
        """
        return cls(x=point.x, y=point.y)

    def separation(self, other):
        r"""Separation to another pixel coordinate.

        This is the two-dimensional cartesian separation :math:`d` with

        .. math::
            d = \sqrt{(x_1 - x_2) ^ 2 + (y_1 - y_2) ^ 2}

        Parameters
        ----------
        other : `PixCoord`
            Other pixel coordinate

        Returns
        -------
        separation : `numpy.array`
            Separation in pixels
        """
        dx = other.x - self.x
        dy = other.y - self.y
        return np.hypot(dx, dy)
