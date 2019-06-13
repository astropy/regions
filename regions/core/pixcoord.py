# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import copy
import numpy as np
from astropy.coordinates import SkyCoord

from .core import _DEFAULT_WCS_ORIGIN, _DEFAULT_WCS_MODE

__all__ = ['PixCoord']


class PixCoord(object):
    """Pixel coordinates.

    This class can represent scalar or array pixel coordinates.

    The data members are either numbers or Numpy arrays
    (not `~astropy.units.Quantity` objects with unit "pixel").

    Given a `astropy.wcs.WCS` object, it can be transformed to and from a
    `~astropy.coordinates.SkyCoord` object.

    Parameters
    ----------
    x : float or array-like
        Pixel coordinate x value
    y : float or array-like
        Pixel coordinate y value

    Examples
    --------

    Usage examples are provided in the :ref:`gs-coord` section of the docs.
    """

    def __init__(self, x, y):

        x, y = np.broadcast_arrays(x, y)

        if x.shape == ():
            self.x, self.y = x.item(), y.item()
        else:
            self.x, self.y = x, y

    def copy(self):
        return self.__class__(
            copy.deepcopy(self.x),
            copy.deepcopy(self.y),
        )

    @staticmethod
    def _validate(val, name, expected='any'):
        """Validate that a given object is an appropriate `PixCoord`.

        This is used for input validation throughout the regions package,
        especially in the `__init__` method of pixel region classes.

        Parameters
        ----------
        val : `PixCoord`
            The object to check
        name : str
            Parameter name (used for error messages)
        expected : {'any', 'scalar', 'not scalar'}
            What kind of PixCoord to check for

        Returns
        -------
        val : `PixCoord`
            The input object (at the moment unmodified, might do fix-ups here later)
        """
        if not isinstance(val, PixCoord):
            raise TypeError('{} must be a PixCoord'.format(name))

        if expected == 'any':
            pass
        elif expected == 'scalar':
            if not val.isscalar:
                raise ValueError('{} must be a scalar PixCoord'.format(name))
        elif expected == 'not scalar':
            if val.isscalar:
                raise ValueError('{} must be a non-scalar PixCoord'.format(name))
        else:
            raise ValueError('Invalid argument for `expected`: {}'.format(expected))

        return val

    @property
    def isscalar(self):
        return np.isscalar(self.x)

    def __repr__(self):
        data = dict(name=self.__class__.__name__, x=self.x, y=self.y)
        fmt = '{name}(x={x}, y={y})'
        return fmt.format(**data)

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

    def __eq__(self, other):
        """
        It checks whether ``other`` is `PixCoord` object and whether their
        abscissa and ordinate values are equal using ``np.assert_allclose``
        with their default tolerance values.
        """
        if isinstance(other, PixCoord):
            return np.allclose([self.x, self.y], [other.x, other.y])
        else:
            return False

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

    @property
    def xy(self):
        """
        A 2-tuple ``(x, y)`` for this coordinate.
        """
        return (self.x, self.y)

    def rotate(self, center, angle):
        """Make a rotated pixel coordinate.

        Rotates counter-clockwise for positive ``angle``.

        Parameters
        ----------
        center : PixCoord
            Rotation center point
        angle : Angle
            Rotation angle

        Returns
        -------
        coord : PixCoord
            Rotated coordinates (an independent copy)
        """
        dx = self.x - center.x
        dy = self.y - center.y
        vec = np.array([dx, dy])

        c, s = np.cos(angle), np.sin(angle)
        rotation_matrix = np.array([[c, -s], [s, c]])

        vec = np.matmul(rotation_matrix, vec)

        return self.__class__(center.x + vec[0], center.y + vec[1])
