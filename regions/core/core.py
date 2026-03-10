# Licensed under a 3-clause BSD style license - see LICENSE.rst
import abc
import copy
import operator

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates.sky_coordinate_parsers import _get_frame_class

from regions._utils.spherical_helpers import bounding_lonlat_poles_processing
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord
from regions.core.registry import RegionsRegistry

__all__ = ['Region', 'PixelRegion', 'SkyRegion',
           'SphericalSkyRegion', 'ComplexSphericalSkyRegion']
__doctest_skip__ = ['Region.serialize', 'Region.write']


def _format_float(val, max_decimals=4, sci_threshold=1e5):
    """
    Format a float for repr/str output.

    Trailing zeros are stripped (a single trailing zero after a decimal
    point is always kept).  Scientific notation is used for values with
    absolute value >= ``sci_threshold`` or <= 1e-5.

    Parameters
    ----------
    val : float
        The float value to format.
    max_decimals : int, optional
        The maximum number of decimal places.
    sci_threshold : float, optional
        The threshold above which scientific notation is used.
    """
    val = float(val)

    if (np.isfinite(val)
            and (abs(val) >= sci_threshold or (0 < abs(val) <= 1e-5))):
        s = f'{val:.{max_decimals}e}'
        mantissa, exponent = s.split('e')
        mantissa = mantissa.rstrip('0')
        if mantissa.endswith('.'):
            mantissa += '0'
        result = f'{mantissa}e{exponent}'
    else:
        result = f'{val:.{max_decimals}f}'
        result = result.rstrip('0')
        if result.endswith('.'):
            result += '0'
    return result


class Region(abc.ABC):
    """
    Base class for all regions.
    """

    _params = ()

    def copy(self, **changes):
        """
        Make an independent (deep) copy.

        Parameters
        ----------
        **changes : dict
            Changes to make to the region parameters.
        """
        fields = list(self._params) + ['meta', 'visual']

        for field in fields:
            if field not in changes:
                changes[field] = copy.deepcopy(getattr(self, field))

        return self.__class__(**changes)

    @staticmethod
    def _format_param_value(val, *, max_decimals=4, max_elements=5):
        """
        Format a region parameter value for use in repr/str output.

        Floats (including pixel coordinates) are formatted using
        `_format_float` with at most ``max_decimals`` decimal
        places and trailing zeros stripped. `~regions.PixCoord`
        arrays and `~astropy.coordinates.SkyCoord` arrays with more
        than ``max_elements`` elements are truncated to the first
        ``max_elements``.

        Parameters
        ----------
        val : object
            The parameter value to format.
        max_decimals : int, optional
            The maximum number of decimal places for floating-point
            values.
        max_elements : int, optional
            The maximum number of array elements to display before
            truncating with ``...``.
        """
        if isinstance(val, SkyCoord):
            if not val.isscalar and len(val) > max_elements:
                result = str(val[:max_elements]).rstrip('>') + ', ...]>'
            else:
                result = str(val)
        elif isinstance(val, PixCoord):
            if val.isscalar:
                if isinstance(val.x, float):
                    x_fmt = _format_float(val.x, max_decimals)
                    y_fmt = _format_float(val.y, max_decimals)
                    result = f'PixCoord(x={x_fmt}, y={y_fmt})'
                else:
                    result = f'PixCoord(x={val.x}, y={val.y})'
            else:
                n = len(val.x)
                xs = val.x[:max_elements] if n > max_elements else val.x
                ys = val.y[:max_elements] if n > max_elements else val.y
                if np.issubdtype(val.x.dtype, np.floating):
                    x_str = ' '.join(
                        _format_float(x, max_decimals) for x in xs)
                    y_str = ' '.join(
                        _format_float(y, max_decimals) for y in ys)
                else:
                    x_str = ' '.join(str(x) for x in xs)
                    y_str = ' '.join(str(y) for y in ys)
                if n > max_elements:
                    x_str += ' ...'
                    y_str += ' ...'
                result = f'PixCoord(x=[{x_str}], y=[{y_str}])'
        elif isinstance(val, u.Quantity):
            result = f'{_format_float(val.value, max_decimals)} {val.unit}'
        elif isinstance(val, float) and not isinstance(val, bool):
            result = _format_float(val, max_decimals)
        else:
            result = str(val)
        return result

    def _cls_info(self):
        """
        Return a list of ``(param_name, formatted_value)`` tuples for
        all region parameters, used by `__repr__` and `__str__`.
        """
        info = []
        if self._params is not None:
            for param in self._params:
                val = getattr(self, param)
                if param == 'text':
                    # Place quotes around text value
                    formatted = repr(val)
                else:
                    formatted = self._format_param_value(val)
                info.append((param, formatted))
        return info

    def __repr__(self):
        prefix = self.__class__.__name__
        cls_info = ', '.join(f'{k}={v}' for k, v in self._cls_info())
        return f'<{prefix}({cls_info})>'

    def __str__(self):
        cls_info = [('Region', self.__class__.__name__)] + self._cls_info()
        return '\n'.join(f'{key}: {val}' for key, val in cls_info)

    def __eq__(self, other):
        """
        Equality operator for Region.

        All Region properties are compared for strict equality except
        for Quantity parameters, which allow for different units if they
        are directly convertible.
        """
        if not isinstance(other, self.__class__):
            return False

        meta_params = ['meta', 'visual']
        self_params = list(self._params) + meta_params
        other_params = list(other._params) + meta_params

        # check that both have identical parameters
        if self_params != other_params:
            return False

        # now check the parameter values
        # Note that Quantity comparisons allow for different units
        # if they directly convertible (e.g., 1. * u.deg == 60. * u.arcmin)
        try:
            for param in self_params:
                # np.any is used for SkyCoord array comparisons
                if np.any(getattr(self, param) != getattr(other, param)):
                    return False
        except TypeError:
            # TypeError is raised from SkyCoord comparison when they do
            # not have equivalent frames. Here return False instead of
            # the TypeError.
            return False

        return True

    def __ne__(self, other):
        """
        Inequality operator for Region.

        Parameters
        ----------
        other : `Region`
            The other region that will be compared.
        """
        return not (self == other)

    @abc.abstractmethod
    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.

        Parameters
        ----------
        other : `Region`
            The other region to use for the intersection.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.

        Parameters
        ----------
        other : `Region`
            The other region to use for the symmetric difference.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.

        Parameters
        ----------
        other : `Region`
            The other region to use for the union.
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

        Parameters
        ----------
        other : `Region`
            The other region to use for the intersection.
        """
        from regions.core.compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other,
                                   operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.

        Parameters
        ----------
        other : `Region`
            The other region to use for the symmetric difference.
        """
        from regions.core.compound import CompoundPixelRegion
        return CompoundPixelRegion(region1=self, region2=other,
                                   operator=operator.xor)

    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.

        Parameters
        ----------
        other : `Region`
            The other region to use for the union.
        """
        from regions.core.compound import CompoundPixelRegion
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

        Notes
        -----
        The conversion between pixel and sky coordinates is an
        approximation. The pixel region shape is mapped to a sky region
        shape using the local pixel scale and angle at the region
        center. Projection effects over the extent of the region are not
        accounted for. The region shape type is always preserved (e.g.,
        a `CirclePixelRegion` converts to a `CircleSkyRegion`).

        For WCS with distortions (e.g., SIP), the local Jacobian matrix
        of the WCS transformation is used to compute directional scale
        factors and angle. For WCS without distortions, a local pixel
        scale and angle are computed using offset-based methods.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_spherical_sky(self, wcs=None, include_boundary_distortions=False,
                         discretize_kwargs=None):
        """
        Convert to an equivalent spherical `~regions.SphericalSkyRegion`
        instance.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance, optional
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates. Required if transforming
            with boundary distortions (if ``include_boundary_distortions`` is True).
            Ignored if boundary distortions not included.

        include_boundary_distortions : bool, optional
            If True, accounts for boundary boundary distortions in spherical to planar
            conversions, by discretizing the boundary and converting the boundary polygon.
            Default is False, which converts to an equivalent idealized shape.

        discretize_kwargs : dict, optional
            Optional keyword arguments to pass to discretize_boundary() method
            if including boundary distortions.

        Returns
        -------
        spherical_sky_region : `~regions.SphericalSkyRegion`
            A spherical sky region, with an equivalent shape (if
            ``include_boundary_distortions`` is False), or a discretized polygon of
            the boundary (if ``include_boundary_distortions`` is True).
        """
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def area(self):
        """
        The exact analytical area of the region shape.
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

        if (mode == 'subpixels'
                and (not isinstance(subpixels, int) or subpixels <= 0)):
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
        artist : `matplotlib.artist.Artist`
            The matplotlib artist (typically a `~matplotlib.patches.Patch`
            object) for the plotted region. The artist can be used, for
            example, when adding a plot legend.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        artist = self.as_artist(origin=origin, **kwargs)
        ax.add_artist(artist)

        return artist


class SkyRegion(Region):
    """
    Base class for all regions defined in celestial coordinates.
    """

    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.

        Parameters
        ----------
        other : `Region`
            The other region to use for the intersection.
        """
        from regions.core.compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other,
                                 operator=operator.and_)

    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.

        Parameters
        ----------
        other : `Region`
            The other region to use for the symmetric difference.
        """
        from regions.core.compound import CompoundSkyRegion
        return CompoundSkyRegion(region1=self, region2=other,
                                 operator=operator.xor)

    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.

        Parameters
        ----------
        other : `Region`
            The other region to use for the union.
        """
        from regions.core.compound import CompoundSkyRegion
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

        Notes
        -----
        The conversion between sky and pixel coordinates is an
        approximation. The sky region shape is mapped to a pixel region
        shape using the local pixel scale and angle at the region
        center. Projection effects over the extent of the region are not
        accounted for. The region shape type is always preserved (e.g.,
        a `CircleSkyRegion` converts to a `CirclePixelRegion`).

        For WCS with distortions (e.g., SIP), the local Jacobian matrix
        of the WCS transformation is used to compute directional scale
        factors and angle. For WCS without distortions, a local pixel
        scale and angle are computed using offset-based methods.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_spherical_sky(self, wcs=None, include_boundary_distortions=False,
                         discretize_kwargs=None):
        """
        Convert to an equivalent spherical `~regions.SphericalSkyRegion`
        instance.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance, optional
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates. Required if transforming
            with boundary distortions (if ``include_boundary_distortions`` is True).
            Ignored if boundary distortions not included.

        include_boundary_distortions : bool, optional
            If True, accounts for boundary boundary distortions in spherical to planar
            conversions, by discretizing the boundary and converting the boundary polygon.
            Default is False, which converts to an equivalent idealized shape.

        discretize_kwargs : dict, optional
            Optional keyword arguments to pass to discretize_boundary() method
            if including boundary distortions.

        Returns
        -------
        spherical_sky_region : `~regions.SphericalSkyRegion`
            A spherical sky region, with an equivalent shape (if
            ``include_boundary_distortions`` is False), or a discretized polygon of
            the boundary (if ``include_boundary_distortions`` is True).
        """
        raise NotImplementedError


class SphericalSkyRegion(Region):
    """
    Base class for all spherical sky regions (compared to the implicitly
    planar SkyRegions).
    """

    def intersection(self, other):
        """
        Return a region representing the intersection of this region
        with ``other``.

        Parameters
        ----------
        other : `~regions.SphericalSkyRegion`
            The other region to use for the intersection.
        """
        from .compound import CompoundSphericalSkyRegion

        return CompoundSphericalSkyRegion(
            region1=self, region2=other, operator=operator.and_
        )

    def symmetric_difference(self, other):
        """
        Return the union of the two regions minus any areas contained in
        the intersection of the two regions.

        Parameters
        ----------
        other : `~regions.SphericalSkyRegion`
            The other region to use for the symmetric difference.
        """
        from .compound import CompoundSphericalSkyRegion

        return CompoundSphericalSkyRegion(
            region1=self, region2=other, operator=operator.xor
        )

    def union(self, other):
        """
        Return a region representing the union of this region with
        ``other``.

        Parameters
        ----------
        other : `~regions.SphericalSkyRegion`
            The other region to use for the union.
        """
        from .compound import CompoundSphericalSkyRegion

        return CompoundSphericalSkyRegion(
            region1=self, region2=other, operator=operator.or_
        )

    @staticmethod
    def _validate_frame(frame):
        # TODO: handle offset origin transformations

        frame = (
            _get_frame_class(frame) if isinstance(frame, str) else frame
        )
        return frame

    @property
    def frame(self):
        """
        Coordinate frame of the region.
        """
        if 'center' in self._params:
            return self.center.frame
        elif 'vertices' in self._params:
            return self.vertices[0].frame
        else:
            raise AttributeError(
                "Either 'center' or 'vertices' must be an attribute/property "
                'of the SphericalSkyRegion.'
            )

    @property
    @abc.abstractmethod
    def bounding_circle(self):
        """
        Bounding circle for the spherical sky region.

        Returns
        -------
        circle_sph_sky_region: `~regions.CircleSphericalSkyRegion`
            A circle spherical sky region object.
        """
        raise NotImplementedError

    def _validate_lonlat_bounds(self, lons_arr, lats_arr, inner_region=None):
        # Check if shape covers either pole & modify lats arr accordingly:
        lons_arr, lats_arr = bounding_lonlat_poles_processing(
            self, lons_arr, lats_arr, inner_region=inner_region
        )
        return lons_arr, lats_arr

    @property
    @abc.abstractmethod
    def bounding_lonlat(self):
        """
        Bounding longitude and latitude of the spherical sky region, in
        the region's frame.

        Returns
        -------
        lons_arr : list of `~astropy.coordinates.Longitude`
            List of lower, upper boundary longitude values.

        lats_arr : list of `~astropy.coordinates.Latitude`
            List of lower, upper boundary latitude values.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def contains(self, coord):
        """
        Check whether a sky coordinate falls inside the spherical sky
        region.

        Parameters
        ----------
        coord : `~astropy.coordinates.SkyCoord`
            The position or positions to check.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def transform_to(self, frame, merge_attributes=True):
        """
        Transform the `~regions.SphericalSkyRegion` instance into
        another instance with a different coordinate reference frame.

        The precise frame transformed to depends on ``merge_attributes``.
        If `False`, the destination frame is used exactly as passed in.
        But this is often not quite what one wants.  E.g., suppose one wants to
        transform an ICRS coordinate that has an obstime attribute to FK4; in
        this case, one likely would want to use this information. Thus, the
        default for ``merge_attributes`` is `True`, in which the precedence is
        as follows: (1) explicitly set (i.e., non-default) values in the
        destination frame; (2) explicitly set values in the source; (3) default
        value in the destination frame.

        Note that in either case, any explicitly set attributes on the source
        |SkyCoord| that are not part of the destination frame's definition are
        kept (stored on the resulting |SkyCoord|), and thus one can round-trip
        (e.g., from FK4 to ICRS to FK4 without losing obstime).

        Parameters
        ----------
        frame : str, or `~astropy.coordinates.BaseCoordinateFrame` class or instance
            The frame to transform this coordinate into.
        merge_attributes : bool, optional
            Whether the default attributes in the destination frame are allowed
            to be overridden by explicitly set attributes in the source
            (see note above; default: `True`).

        Returns
        -------
        sph_sky_region : `~regions.SphericalSkyRegion`
            A new spherical sky region represented in the ``frame`` frame.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def discretize_boundary(self, n_points=10):
        """
        Discretize the boundary into a
        `~regions.PolygonSphericalSkyRegion`, as an approximation where
        all sides follow great circles.

        Parameters
        ----------
        n_points : int, optional
            Number of points along the region's boundary.

        Returns
        -------
        poly_sky_region: `~regions.PolygonSphericalSkyRegion`
            Spherical sky polygon object.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_sky(
        self, wcs=None, include_boundary_distortions=False, discretize_kwargs=None
    ):
        """
        Convert to a planar `~regions.SkyRegion` instance.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance, optional
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates. Required if transforming
            with boundary distortions (if ``include_boundary_distortions`` is True).
            Ignored if boundary distortions not included.

        include_boundary_distortions : bool, optional
            If True, accounts for boundary boundary distortions in spherical to planar
            conversions, by discretizing the boundary and converting the boundary polygon.
            Default is False, which converts to an equivalent idealized shape.

        discretize_kwargs : dict, optional
            Optional keyword arguments to pass to discretize_boundary() method
            if including boundary distortions.

        Returns
        -------
        sky_region : `~regions.SkyRegion`
            A planar sky region, with an equivalent shape (if
            ``include_boundary_distortions`` is False), or a discretized polygon of
            the boundary (if ``include_boundary_distortions`` is True).
        """
        raise NotImplementedError

    @abc.abstractmethod
    def to_pixel(
        self, wcs=None, include_boundary_distortions=False, discretize_kwargs=None
    ):
        """
        Convert to a planar `~regions.PixelRegion` instance.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS` instance, optional
            The world coordinate system transformation to use to convert
            between sky and pixel coordinates. Required if transforming
            with boundary distortions (if ``include_boundary_distortions`` is True).
            Ignored if boundary distortions not included.

        include_boundary_distortions : bool, optional
            If True, accounts for boundary boundary distortions in spherical to planar
            conversions, by discretizing the boundary and converting the boundary polygon.
            Default is False, which converts to an equivalent idealized shape.

        discretize_kwargs : dict, optional
            Optional keyword arguments to pass to discretize_boundary() method
            if including boundary distortions.

        Returns
        -------
        pixel_region : `~regions.PixelRegion`
            A pixel region, with an equivalent shape (if
            ``include_boundary_distortions`` is False), or a discretized polygon of
            the boundary (if ``include_boundary_distortions`` is True).
        """
        raise NotImplementedError

    def write(self, filename, format=None, overwrite=False, **kwargs):
        # Note explicitly that this is not yet implemented.
        raise NotImplementedError

    def serialize(self, format=None, **kwargs):
        # Note explicitly that this is not yet implemented.
        raise NotImplementedError


class ComplexSphericalSkyRegion(SphericalSkyRegion):
    """
    Base class for complex cases, where the definitional parameters do
    not / cannot transform between coordinate frames (including
    `~regions.RangeSphericalSkyRegion`).
    """

    # Because the parameters don't transform,
    # these object cannot be defined simply as "compound regions"

    # So some repr and str methods are changed.

    # (For now, this is just RangeSphericalSkyRegion,
    # but breaking out these methods/properties to streamline).

    # Instead, a variable set _boundaries is used, containing the
    # subclass boundary properties. These boundary properties are, in order of preference,
    # (a) accessed from `_[boundname]`, for transformed instances of subclass shapes that
    #     don't simply map in parameterization from frame to frame (eg, Range)
    # (b) derived on-the-fly

    _boundaries = ()

    def __repr__(self):
        prefix = f'{self.__class__.__name__}'
        cls_info = []

        _do_params_info = False
        if (self._params is not None) and ((len(self._params) > 1) | (self._params[0] != 'frame')):
            # Only do param-based repr if params other than "frame" are set
            # If only param is "frame", make repr with boundary info
            _do_params_info = True

        if _do_params_info:
            for param in self._params:
                if param == 'text':
                    # place quotes around text value
                    keyval = f'{param}={getattr(self, param)!r}'
                elif param == 'frame':
                    attr = getattr(self, param)
                    keyval = f'{param}={getattr(attr, "name", repr(attr))}'
                else:
                    keyval = f'{param}={getattr(self, param)}'
                cls_info.append(keyval)

            cls_info = ',\n'.join(cls_info)
            return f'<{prefix}(\n{cls_info}\n)>'
        else:
            # First check if "frame" in self._params:
            if (self._params is not None) and self._params[0] == 'frame':
                param = 'frame'
                attr = getattr(self, param)
                keyval = f'{param}={getattr(attr, "name", repr(attr))}'
                cls_info.append(keyval)

            # If "params" is None, eg for a "transformed"
            # complex shape which is no longer described by the "high-level"
            # region parameters, instead give information on the regions
            # contained within _boundaries info:
            for bound in self._boundaries:
                keyval = f'{bound}={repr(getattr(self, bound))}'
                cls_info.append(keyval)

            cls_info = ',\n'.join(cls_info)
            return f'<{prefix}(\n{cls_info}\n)>'

    def __str__(self):
        cls_info = [('Region', self.__class__.__name__)]

        _do_params_info = False
        if (self._params is not None) and ((len(self._params) > 1) | (self._params[0] != 'frame')):
            # Only do param-based str if params other than "frame" are set.
            # If only param is "frame", make str with boundary info
            _do_params_info = True

        if _do_params_info:
            for param in self._params:
                if param in ['text']:
                    keyval = (param, repr(getattr(self, param)))

                elif param == 'frame':
                    attr = getattr(self, param)
                    keyval = (param, getattr(attr, 'name', repr(attr)))
                else:
                    keyval = (param, getattr(self, param))
                cls_info.append(keyval)
            return '\n'.join([f'{key}: {val}' for key, val in cls_info])
        else:
            # First check if "frame" in self._params:
            if (self._params is not None) and (self._params[0] == 'frame'):
                param = 'frame'
                attr = getattr(self, param)
                keyval = (param, getattr(attr, 'name', repr(attr)))
                cls_info.append(keyval)

            # If "params" is None, eg for a transformed complex shape,
            # instead give information on the regions contained
            # within _boundaries info:
            for bound in self._boundaries:
                keyval = (bound, repr(getattr(self, bound)))
                cls_info.append(keyval)

            return '\n'.join([f'{key}: {val}' for key, val in cls_info])

    @property
    @abc.abstractmethod
    def _compound_region(self):
        """
        Compound region containing composite boundaries for this complex
        region shape.
        """

    def copy(self, **changes):
        # Boundaries: only copy stashed internal attributes "_{bound name}",
        # as otherwise these are derived on-the-fly.
        boundaries_interalattr = [f"_{bn}" for bn in list(self._boundaries)]
        fields = boundaries_interalattr + [
            '_frame', '_vertices',
            '_is_original_frame', '_params',
        ] + ['meta', 'visual']
        if self._params is not None:
            fields += list(self._params)

        for field in fields:
            if (field not in changes) & hasattr(self, field):
                changes[field] = copy.deepcopy(getattr(self, field))

        return self.__class__(**changes)

    @property
    def frame(self):
        frame = getattr(self, '_frame', None)
        if frame is None:
            return self._compound_region.frame
        return frame
