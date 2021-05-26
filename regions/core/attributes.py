# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The module provides several custom descriptor classes for attribute
validation of region classes.
"""

import abc

from astropy.coordinates import SkyCoord
from astropy.units import Quantity
import numpy as np

from .core import PixelRegion, SkyRegion
from .pixcoord import PixCoord

__all__ = []


class RegionAttr(abc.ABC):
    """
    Base descriptor class for region attribute validation.
    """

    def __init__(self, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return instance.__dict__[self.name]

    def __set__(self, instance, value):
        self._validate(value)
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        del instance.__dict__[self.name]

    @abc.abstractmethod
    def _validate(self, value):
        """
        Validate the attribute value.

        An exception is raised if the value is invalid.
        """
        pass


class ScalarPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a scalar
    `~regions.PixCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and value.isscalar):
            raise ValueError(f'The {self.name} must be a scalar PixCoord '
                             'object')


class OneDPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a
    one-dimensional `regions.PixCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and not value.isscalar
                and value.x.ndim == 1):
            raise ValueError(f'The {self.name} must be a 1D PixCoord object')


class ScalarLength(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a scalar
    python/numpy number.
    """

    def _validate(self, value):
        if not np.isscalar(value):
            raise ValueError(
                f'The {self.name} must be a scalar numpy/python number')


class ScalarSky(RegionAttr):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a scalar
    `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.isscalar):
            raise ValueError(f'The {self.name} must be a scalar SkyCoord '
                             'object')


class OneDSky(RegionAttr):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a
    one-dimensional `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.ndim == 1):
            raise ValueError(f'The {self.name} must be a 1D SkyCoord object')


class QuantityLength(RegionAttr):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a scalar
    `~astropy.units.Quantity` object.
    """

    def _validate(self, value):
        if not (isinstance(value, Quantity) and value.isscalar):
            raise ValueError(f'The {self.name} must be a scalar astropy '
                             'Quantity object')


class CompoundRegionPix(RegionAttr):
    """
    Descriptor class for `~regions.CompoundPixelRegion`, which takes a
    `~regions.PixelRegion` object.
    """

    def _validate(self, value):
        if not isinstance(value, PixelRegion):
            raise ValueError(f'The {self.name} must be a PixelRegion object')


class CompoundRegionSky(RegionAttr):
    """
    Descriptor class for `~regions.CompoundSkyRegion`, which takes a
    `~regions.SkyRegion` object.
    """

    def _validate(self, value):
        if not isinstance(value, SkyRegion):
            raise ValueError(f'The {self.name} must be a SkyRegion object')
