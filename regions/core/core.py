# Licensed under a 3-clause BSD style license - see LICENSE.rst
import abc
import copy
import operator

from .metadata import RegionMeta, RegionVisual
from .pixcoord import PixCoord
from .registry import RegionsRegistry

__all__ = ['Region', 'PixelRegion', 'SkyRegion']
__doctest_skip__ = ['Region.serialize', 'Region.write']


class Region(abc.ABC):
    """
    Base class for all regions.
    """

    _params = ()

    def copy(self, **changes):
        """
        Make an independent (deep) copy.
        """
        fields = list(self._params) + ['meta', 'visual']

        for field in fields:
            if field not in changes:
                changes[field] = copy.deepcopy(getattr(self, field))

        return self.__class__(**changes)

    def __repr__(self):
        prefix = f'{self.__class__.__name__}'
        cls_info = []
        if self._params is not None:
            for param in self._params:
                cls_info.append(f'{param}={getattr(self, param)}')
        cls_info = ', '.join(cls_info)
        return f'<{prefix}({cls_info})>'

    def __str__(self):
        cls_info = [('Region', self.__class__.__name__)]
        if self._params is not None:
            for param in self._params:
                cls_info.append((param, getattr(self, param)))
        return '\n'.join([f'{key}: {val}' for key, val in cls_info])

    @abc.abstractmethod
    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.
        """
        raise NotImplementedError

    def __and__(self, other):
        return self.intersection(other)

    def __or__(self, other):
        return self.union(other)

    def __xor__(self, other):
        return self.symmetric_difference(other)

    @classmethod
    def get_formats(cls):
        """
        Get the registered I/O formats as a Table.
        """
        return RegionsRegistry.get_formats(Region)

    def write(self, filename, format=None, overwrite=False, **kwargs):
        """
        Write the region to a region file in the specified format.

        This method allows writing a file in many supported data
        formats, e.g.,::

            >>> reg.write('new_regions.reg', format='ds9')
            >>> reg.write('new_regions.crtf', format='crtf')
            >>> reg.write('new_regions.fits', format='fits')

        A list of the available formats for `~regions.Region` is
        available using::

            >>> from regions import Region
            >>> Region.get_formats()

        Parameters
        ----------
        filename : str
            The filename or URL of the file to write.

        format : str, optional
            The file format specifier.

        overwrite : bool, optional
            If True, overwrite the output file if it exists. Raises an
            `OSError` if False and the output file exists. Default is
            False.

        **kwargs : dict, optional
            Keyword arguments passed to the data writer.
        """
        return RegionsRegistry.write([self], filename,
                                     Region, format=format,
                                     overwrite=overwrite, **kwargs)

    def serialize(self, format=None, **kwargs):
        """
        Serialize the region to a region string or table.

        This method allows serializing regions in many supported data
        formats, e.g.,::

            >>> reg1_str = reg.serialize(format='ds9')
            >>> reg2_str = reg.serialize(format='crtf')
            >>> reg3_tbl = reg.serialize(format='fits')

        A list of the available formats for `~regions.Region` is
        available using::

            >>> from regions import Region
            >>> Region.get_formats()

        Parameters
        ----------
        format : str, optional
            The file format specifier.

        **kwargs : dict, optional
            Keyword arguments passed to the data serializer.
        """
        return RegionsRegistry.serialize([self], Region, format=format,
                                         **kwargs)


class PixelRegion(Region):
    """
    Base class for all regions defined in pixel coordinates.
    """

    meta = RegionMeta()
    visual = RegionVisual()

    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other,
                                   operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other,
                                   operator=operator.xor)

    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.
        """
        from .compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other,
                                   operator=operator.or_)

    @abc.abstractmethod
    def contains(self, pixcoord):
        """
        Check whether a position or positions fall inside the region.

        Parameters
        ----------
        pixcoord : `~regions.PixCoord`
            The position or positions to check.
        """
        raise NotImplementedError

    def __contains__(self, coord):
        if not coord.isscalar:
            raise ValueError(f'coord must be scalar. coord={coord}')
        return self.contains(coord)

    @abc.abstractmethod
    def to_sky(self, wcs):
        """
        Return a region defined in sky coordinates.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The world coordinate system transformation to use to convert
            from pixels to sky coordinates.

        Returns
        -------
        sky_region : `~regions.SkyRegion`
            The sky region.
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def area(self):
        """
        The area of the region.
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def bounding_box(self):
        """
        The minimal bounding box (in integer pixel coordinates) that
        contains the region.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_mask(self, mode='center', subpixels=5):
        """
        Return a mask for the region.

        Parameters
        ----------
        mode : {'center', 'exact', 'subpixels'}, optional
            The method used to determine the overlap of the region on
            the pixel grid. Not all options are available for all region
            types. Note that the more precise methods are generally
            slower. The following methods are available:

                * ``'center'``:
                  A pixel is considered to be entirely in or out of the
                  region depending on whether its center is in or out of
                  the region. The returned mask will contain values only
                  of 0 (out) and 1 (in).

                * ``'exact'`` (default):
                  The exact fractional overlap of the region and each
                  pixel is calculated. The returned mask will contain
                  values between 0 and 1.

                * ``'subpixel'``:
                  A pixel is divided into subpixels (see the
                  ``subpixels`` keyword), each of which are considered
                  to be entirely in or out of the region depending
                  on whether its center is in or out of the region.
                  If ``subpixels=1``, this method is equivalent to
                  ``'center'``. The returned mask will contain values
                  between 0 and 1.

        subpixels : int, optional
            For the ``'subpixel'`` mode, resample pixels by this factor
            in each dimension. That is, each pixel is divided into
            ``subpixels ** 2`` subpixels.

        Returns
        -------
        mask : `~regions.RegionMask`
            A mask for the region.
        """
        raise NotImplementedError

    @staticmethod
    def _validate_mode(mode, subpixels):
        valid_modes = ('center', 'exact', 'subpixels')
        if mode not in valid_modes:
            raise ValueError(f'Invalid mask mode: {mode} (should be one '
                             f'of {valid_modes}')

        if mode == 'subpixels':
            if not isinstance(subpixels, int) or subpixels <= 0:
                raise ValueError(f'Invalid subpixels value: {subpixels} '
                                 '(should be a strictly positive integer)')

    @abc.abstractmethod
    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Convert to matplotlib patch object for this region.

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Patch`.

        Returns
        -------
        patch : `~matplotlib.patches.Patch`
            A matplotlib patch.
        """
        raise NotImplementedError

    def _define_mpl_kwargs(self, artist='Patch'):
        """
        Define a dictionary of matplotlib keywords for the input
        ``artist`` from the region's ``visual`` properties.

        Unset ``visual`` properties are set to default values based on
        the DS9 convention.
        """
        if artist not in ('Patch', 'Line2D', 'Text'):
            raise ValueError(f'artist "{artist}" is not supported')

        kwargs = {}
        kwargs['color'] = self.visual.get('color', 'green')
        kwargs['label'] = self.meta.get('label', '')

        if artist == 'Text':
            kwargs['family'] = self.visual.get('font', 'helvetica')
            kwargs['style'] = self.visual.get('fontsyle', 'normal')
            kwargs['weight'] = self.visual.get('fontweight', 'roman')
            kwargs['size'] = self.visual.get('fontsize', '12')
            kwargs['rotation'] = self.visual.get('textangle', '0')
        else:
            kwargs['linewidth'] = self.visual.get('linewidth', 1)
            kwargs['linestyle'] = self.visual.get('linstyle', 'solid')

            if artist == 'Line2D':
                kwargs['marker'] = self.visual.get('symbol', 'o')
                kwargs['markersize'] = self.visual.get('symsize', 11)
                kwargs['markeredgecolor'] = kwargs['color']
                kwargs['markeredgewidth'] = self.visual.get('width', 1)
                kwargs['fillstyle'] = self.visual.get('fill', 'none')

            if artist == 'Patch':
                kwargs['edgecolor'] = kwargs.pop('color')
                kwargs['fill'] = self.visual.get('fill', False)

        return kwargs

    def plot(self, origin=(0, 0), ax=None, **kwargs):
        """
        Plot the region on a matplotlib `~matplotlib.axes.Axes`
        instance.

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        ax : `~matplotlib.axes.Axes` or `None`, optional
            The matplotlib axes on which to plot.  If `None`, then the
            current `~matplotlib.axes.Axes` instance is used.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Patch`.

        Returns
        -------
        ax : `~matplotlib.axes.Axes`
            The axes on which the patch was added.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        artist = self.as_artist(origin=origin, **kwargs)
        ax.add_artist(artist)

        return ax


class SkyRegion(Region):
    """
    Base class for all regions defined in celestial coordinates.
    """

    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other,
                                 operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other,
                                 operator=operator.xor)

    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.
        """
        from .compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other,
                                 operator=operator.or_)

    def contains(self, skycoord, wcs):
        """
        Check whether a sky coordinate falls inside the region.

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The position or positions to check.
        wcs : `~astropy.wcs.WCS`
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates.
        """
        pixel_region = self.to_pixel(wcs)
        pixcoord = PixCoord.from_sky(skycoord, wcs)
        return pixel_region.contains(pixcoord)

    @abc.abstractmethod
    def to_pixel(self, wcs):
        """
        Return the equivalent region defined in pixel coordinates.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates.

        Returns
        -------
        pixel_region : `~regions.PixelRegion`
            A pixel region.
        """
        raise NotImplementedError
