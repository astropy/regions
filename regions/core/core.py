# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import abc
import operator

from astropy.extern import six


__all__ = ['Region', 'PixelRegion', 'SkyRegion']


"""
Here we define global variables for the default `origin` and `mode` used
for WCS transformations throughout the `regions` package.

Their purpose is to simplify achieving uniformity across the codebase.
They are mainly used as default arguments for methods that do WCS
transformations.

They are private (with an underscore), not part of the public API,
users should not touch them.
"""
_DEFAULT_WCS_ORIGIN = 0
_DEFAULT_WCS_MODE = 'all'

VALID_MASK_MODES = set(['center', 'exact', 'subpixels'])


@six.add_metaclass(abc.ABCMeta)
class Region(object):
    """
    Base class for all regions.
    """

    def __repr__(self):
        if hasattr(self, 'center'):
            params = [repr(self.center)]
        else:
            params = []
        if self._repr_params is not None:
            for key, val in self._repr_params:
                params.append('{0}={1}'.format(key, val))
        params = ', '.join(params)

        return '<{0}({1})>'.format(self.__class__.__name__, params)

    def __str__(self):
        cls_info = [('Region', self.__class__.__name__)]
        if hasattr(self, 'center'):
            cls_info.append(('center', self.center))
        if self._repr_params is not None:
            cls_info += self._repr_params
        fmt = ['{0}: {1}'.format(key, val) for key, val in cls_info]

        return '\n'.join(fmt)

    @abc.abstractmethod
    def intersection(self, other):
        """
        Returns a region representing the intersection of this region with
        ``other``.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def symmetric_difference(self, other):
        """
        Returns the union of the two regions minus any areas contained in the
        intersection of the two regions.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def union(self, other):
        """
        Returns a region representing the union of this region with ``other``.
        """
        raise NotImplementedError

    def __and__(self, other):
        return self.intersection(other)

    def __or__(self, other):
        return self.union(other)

    def __xor__(self, other):
        return self.symmetric_difference(other)

    def __contains__(self, coord):
        if not coord.isscalar:
            raise ValueError('coord must be scalar. coord={}'.format(coord))
        return self.contains(coord)


@six.add_metaclass(abc.ABCMeta)
class PixelRegion(Region):
    """
    Base class for all regions defined in pixel coordinates
    """

    def intersection(self, other):
        """
        Returns a region representing the intersection of this region with
        ``other``.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(self, other, operator.and_)

    def symmetric_difference(self, other):
        """
        Returns the union of the two regions minus any areas contained in the
        intersection of the two regions.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(self, other, operator.xor)

    def union(self, other):
        """
        Returns a region representing the union of this region with ``other``.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(self, other, operator.or_)

    @abc.abstractmethod
    def contains(self, pixcoord):
        """
        Checks whether a position or positions fall inside the region.

        Parameters
        ----------
        pixcoord : `~regions.PixCoord`
            The position or positions to check, as a tuple of scalars or
            arrays. In future this could also be a `PixCoord` instance.
        """

    @abc.abstractmethod
    def to_sky(self, wcs, mode='local', tolerance=None):
        """
        Returns a region defined in sky coordinates.

        Parameters
        ----------

        wcs : `~astropy.wcs.WCS` instance
            The world coordinate system transformation to assume

        mode : str
            Converting to sky coordinates can be done with various degrees of
            approximation, which can be set with this option. Possible values
            are:

            * `'local'`: assume that the field of view is small and that
              pixels are square, so that e.g. a circle in sky coordinates
              would be a circle in pixel coordinates. This is the fastest and
              most commonly used for e.g. photometry.

            * `'affine'`: approximate any deviations from the 'local'
              assumption by an affine transformation, e.g. a circle would
              become a rotated ellipse.

            * `'full'`: return an arbitrarily complex polygon in sky
              coordinates that represents the full level of distortion due to
              the conversion from pixel to world coordinates. The degree of
              exactness can be controlled by the ``tolerance`` argument.

        tolerance : `~astropy.units.Quantity`
            The tolerance for the ``'full'`` mode described above.
        """

    @abc.abstractproperty
    def bounding_box(self):
        """
        The minimal bounding box (in integer pixel coordinates) that contains
        the region.
        """
        pass

    @abc.abstractmethod
    def to_mask(self, mode='center', subpixels=5):
        """
        Returns a mask for the aperture.

        Parameters
        ----------
        mode : { 'center' | 'exact' | 'subpixels'}
            The following modes are available:
                * ``'center'``: returns 1 for pixels where the center is in
                  the region, and 0 otherwise.
                * ``'exact'``: returns a value between 0 and 1 giving the
                  fractional level of overlap of the pixel with the region.
                * ``'subpixels'``: A pixel is divided into subpixels and
                  the center of each subpixel is tested (a subpixel is
                  either completely in or out of the region).  Returns a
                  value between 0 and 1 giving the fractional level of
                  overlap of the subpixels with the region.  With
                  ``subpixels`` set to 1, this method is equivalent to
                  ``'center'``.
        subpixels : int, optional
            For the ``'subpixel'`` mode, resample pixels by this factor
            in each dimension. That is, each pixel is divided into
            ``subpixels ** 2`` subpixels.

        Returns
        -------
        mask : `~regions.Mask`
            A region mask object.
        """

    @staticmethod
    def _validate_mode(mode, subpixels):
        if mode not in VALID_MASK_MODES:
            raise ValueError("Invalid mask mode: {0} (should be one "
                             "of {1})".format(mode, '/'.join(VALID_MASK_MODES)))
        if mode == 'subpixels':
            if not isinstance(subpixels, int) or subpixels <= 0:
                raise ValueError("Invalid subpixels value: {0} (should be"
                                 " a strictly positive integer)".format(subpixels))

    @abc.abstractmethod
    def to_shapely(self):
        """
        Convert this region to a Shapely object.
        """

    @abc.abstractmethod
    def as_patch(self, **kwargs):
        """Convert to mpl patch

        Returns
        -------
        patch : `~matplotlib.patches.Patch`
            Matplotlib patch
        """

    def plot(self, ax=None, **kwargs):
        """
        Calls as_patch method forwarding all kwargs and adds patch
        to given axis.

        Parameters
        ----------
        ax : `~matplotlib.axes`, optional
            Axis
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        patch = self.as_patch(**kwargs)
        ax.add_patch(patch)

        return ax


@six.add_metaclass(abc.ABCMeta)
class SkyRegion(Region):
    """
    Base class for all regions defined in celestial coordinates
    """

    def intersection(self, other):
        """
        Returns a region representing the intersection of this region with
        ``other``.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(self, other, operator.and_)

    def symmetric_difference(self, other):
        """
        Returns the union of the two regions minus any areas contained in the
        intersection of the two regions.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(self, other, operator.xor)

    def union(self, other):
        """
        Returns a region representing the union of this region with ``other``.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(self, other, operator.or_)

    @abc.abstractmethod
    def contains(self, skycoord):
        """
        Checks whether a position or positions fall inside the region.

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The position or positions to check
        """

    @abc.abstractmethod
    def to_pixel(self, wcs, mode='local', tolerance=None):
        """
        Returns a region defined in pixel coordinates.

        Parameters
        ----------

        wcs : `~astropy.wcs.WCS` instance
            The world coordinate system transformation to assume

        mode : str
            Converting to pixel coordinates can be done with various degrees
            of approximation, which can be set with this option. Possible
            values are:

            * `'local'`: assume that the field of view is small and that
              pixels are square, so that e.g. a circle in sky coordinates
              would be a circle in pixel coordinates. This is the fastest and
              most commonly used for e.g. photometry.

            * `'affine'`: approximate any deviations from the 'local'
              assumption by an affine transformation, e.g. a circle would
              become a rotated ellipse.

            * `'full'`: return an arbitrarily complex polygon in pixel
              coordinates that represents the full level of distortion due to
              the conversion from world to pixel coordinates. The degree of
              exactness can be controlled by the ``tolerance`` argument.

        tolerance : `~astropy.units.Quantity`
            The tolerance for the ``'full'`` mode described above.
        """
