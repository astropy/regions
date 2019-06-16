# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import abc
import six
import copy
import operator
import inspect


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

VALID_MASK_MODES = {'center', 'exact', 'subpixels'}


class MetaRegion(abc.ABCMeta):
    """
    This is a sub metaclass of `abc.ABCMeta` that makes methods of a class
    automatically have their docstrings filled in from the methods they override
    in the base class.

    If the class uses multiple inheritance, the docstring will be
    chosen from the first class in the bases list, in the same way as
    methods are normally resolved in Python.  If this results in
    selecting the wrong docstring, the docstring will need to be
    explicitly included on the method.

    """

    def __init__(cls, name, bases, dct):
        def is_public_member(key):
            return (
                (key.startswith('__') and key.endswith('__')
                 and len(key) > 4) or
                not key.startswith('_'))

        for key, val in six.iteritems(dct):
            if (inspect.isfunction(val) and is_public_member(key)
                    and val.__doc__ is None):
                for base in cls.__mro__[1:]:
                    super_method = getattr(base, key, None)
                    if super_method is not None:
                        val.__doc__ = super_method.__doc__
                        break

        super(MetaRegion, cls).__init__(name, bases, dct)


@six.add_metaclass(MetaRegion)
class Region(object):
    """
    Base class for all regions.
    """

    def copy(self, **changes):
        """Make an independent (deep) copy."""
        fields = []
        if hasattr(self, "center"):
            fields.append("center")
        if hasattr(self, "vertices"):
            fields.append("vertices")
        fields.extend(self._repr_params)
        fields.extend(["meta", "visual"])

        for field in fields:
            if field not in changes:
                changes[field] = copy.deepcopy(getattr(self, field))

        return self.__class__(**changes)

    def __repr__(self):
        if hasattr(self, 'center'):
            params = [repr(self.center)]
        else:
            params = []
        if self._repr_params is not None:
            for key in self._repr_params:
                params.append('{0}={1}'.format(key.replace("_", " "),
                                               getattr(self, key)))
        params = ', '.join(params)

        return '<{0}({1})>'.format(self.__class__.__name__, params)

    def __str__(self):
        cls_info = [('Region', self.__class__.__name__)]
        if hasattr(self, 'center'):
            cls_info.append(('center', self.center))
        if self._repr_params is not None:
            params_value = [(x.replace("_", " "), getattr(self, x))
                            for x in self._repr_params]
            cls_info += params_value
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
        return CompoundPixelRegion(region1=self, region2=other, operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Returns the union of the two regions minus any areas contained in the
        intersection of the two regions.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other, operator=operator.xor)

    def union(self, other):
        """
        Returns a region representing the union of this region with ``other``.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other, operator=operator.or_)

    @abc.abstractmethod
    def contains(self, pixcoord):
        """
        Checks whether a position or positions fall inside the region.

        Parameters
        ----------
        pixcoord : `~regions.PixCoord`
            The position or positions to check.
        """
        raise NotImplementedError

    def __contains__(self, coord):
        if not coord.isscalar:
            raise ValueError('coord must be scalar. coord={}'.format(coord))
        return self.contains(coord)

    @abc.abstractmethod
    def to_sky(self, wcs):
        """
        Returns a region defined in sky coordinates.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance
            The world coordinate system transformation to assume

        Returns
        -------
        sky_region : `~regions.SkyRegion` object.
        """
        raise NotImplementedError

    @abc.abstractproperty
    def bounding_box(self):
        """
        The minimal bounding box (in integer pixel coordinates) that contains
        the region.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_mask(self, mode='center', subpixels=5):
        """
        Returns a mask for the aperture.

        Parameters
        ----------
        mode : { 'center' | 'exact' | 'subpixels'}, optional
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
        raise NotImplementedError

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
    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Convert to mpl patch

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).

        kwargs : `dict`
            keywords that a `~matplotlib.patches.Patch` accepts

        Returns
        -------
        patch : `~matplotlib.patches.Patch`
            Matplotlib patch
        """
        raise NotImplementedError

    def mpl_properties_default(self, shape='patch'):
        """
        This sets the default values of the visual attributes as specified
        under DS9 convention.
        """

        kwargs = dict()
        kwargs['color'] = self.visual.get('color', 'green')
        kwargs['label'] = self.meta.get('label', "")

        if shape == 'text':
            kwargs['family'] = self.visual.get('font', 'helvetica')
            kwargs['rotation'] = self.visual.get('textangle', '0')
            kwargs['size'] = self.visual.get('fontsize', '12')
            kwargs['style'] = self.visual.get('fontsyle', 'normal')
            kwargs['weight'] = self.visual.get('fontweight', 'roman')

        else:
            if shape == 'Line2D':
                kwargs['marker'] = self.visual.get('symbol', 'o')
                kwargs['markersize'] = self.visual.get('symsize', 11)
                kwargs['markeredgecolor'] = kwargs['color']
                kwargs['markeredgewidth'] = self.visual.get('width', 1)
            if shape == 'patch':
                kwargs['edgecolor'] = kwargs.pop('color')
                kwargs['fill'] = self.visual.get('fill', False)

            kwargs['linewidth'] = self.visual.get('linewidth', 1)
            kwargs['linestyle'] = self.visual.get('linstyle', 'solid')

        return kwargs

    def plot(self, origin=(0, 0), ax=None, **kwargs):
        """
        Calls ``as_artist`` method forwarding all kwargs and adds patch
        to given axis.

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed image.
            Default is (0, 0).
        ax : `~matplotlib.axes.Axes`, optional
            Axis
        kwargs : `dict`
            keywords that a `~matplotlib.patches.Patch` accepts

        Returns
        -------
        ax : `~matplotlib.axes.Axes`
            Axes on which the patch is added.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        artist = self.as_artist(origin=origin, **kwargs)
        ax.add_artist(artist)

        return ax


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
        return CompoundSkyRegion(region1=self, region2=other, operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Returns the union of the two regions minus any areas contained in the
        intersection of the two regions.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other, operator=operator.xor)

    def union(self, other):
        """
        Returns a region representing the union of this region with ``other``.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other, operator=operator.or_)

    def contains(self, skycoord, wcs):
        """
        Check whether a sky coordinate falls inside the region

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The position or positions to check
        wcs : `~astropy.wcs.WCS` instance
            The world coordinate system transformation to assume
        """
        from .pixcoord import PixCoord
        pixel_region = self.to_pixel(wcs)
        pixcoord = PixCoord.from_sky(skycoord, wcs)
        return pixel_region.contains(pixcoord)

    @abc.abstractmethod
    def to_pixel(self, wcs):
        """
        Returns the equivalent region defined in pixel coordinates.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance
            The world coordinate system transformation to assume

        Returns
        -------
        pixel_region : `~regions.PixelRegion` object.
        """
        raise NotImplementedError
