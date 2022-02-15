# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The module provides several custom descriptor classes for attribute
validation of region classes.
"""

import abc

from astropy.coordinates import SkyCoord
from astropy.units import Quantity
import numpy as np

from .pixcoord import PixCoord

__all__ = []


class RegionAttribute(abc.ABC):
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


class ScalarPix(RegionAttribute):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a scalar
    `~regions.PixCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and value.isscalar):
            raise ValueError(f'The {self.name} must be a scalar PixCoord '
                             'object')


class OneDPix(RegionAttribute):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a
    one-dimensional `regions.PixCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and not value.isscalar
                and value.x.ndim == 1):
            raise ValueError(f'The {self.name} must be a 1D PixCoord object')


class PositiveScalar(RegionAttribute):
    """
    Descriptor class to check that value is a strictly positive (> 0)
    scalar.
    """

    def _validate(self, value):
        if not np.isscalar(value) or value <= 0:
            raise ValueError(f'{self.name} must be a positive scalar')


class ScalarSky(RegionAttribute):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a scalar
    `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.isscalar):
            raise ValueError(f'The {self.name} must be a scalar SkyCoord '
                             'object')


class OneDSky(RegionAttribute):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a
    one-dimensional `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.ndim == 1):
            raise ValueError(f'The {self.name} must be a 1D SkyCoord object')


class ScalarAngle(RegionAttribute):
    """
    Descriptor class to check that value is a scalar angle, either an
    astropy Angle or Quantity with angular units.
    """

    def _validate(self, value):
        if isinstance(value, Quantity):
            if not value.isscalar:
                raise ValueError(f'{self.name} must be a scalar')

            if not value.unit.physical_type == 'angle':
                raise ValueError(f'{self.name} must have angular units')
        else:
            raise ValueError(f'{self.name} must be a scalar angle')


class PositiveScalarAngle(RegionAttribute):
    """
    Descriptor class to check that value is a strictly positive scalar
    angle, either an astropy Angle or Quantity with angular units.
    """

    def _validate(self, value):
        if isinstance(value, Quantity):
            if not value.isscalar:
                raise ValueError(f'{self.name} must be a scalar')

            if not value.unit.physical_type == 'angle':
                raise ValueError(f'{self.name} must have angular units')

            if not value > 0:
                raise ValueError(f'{self.name} must be strictly positive')
        else:
            raise ValueError(f'{self.name} must be a positive scalar angle')


class RegionType(RegionAttribute):
    """
    Descriptor class for compound pixel and sky regions.
    """

    def __init__(self, name, regionclass):
        self.name = name
        self.regionclass = regionclass

    def _validate(self, value):
        if not isinstance(value, self.regionclass):
            raise ValueError(f'The {self.name} must be a '
                             f'{self.regionclass.__name__} object')
