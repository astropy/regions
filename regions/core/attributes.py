# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The module provides several custom descriptor classes for attribute validation of region classes.  It also contains RegionMeta and RegionVisual classes to handle meta data of regions.
"""

import abc

from astropy.coordinates import SkyCoord
from astropy.units import Quantity

import numpy as np

from .pixcoord import PixCoord
from .core import PixelRegion, SkyRegion

__all__ = ['Meta', 'RegionMeta', 'RegionVisual']


class RegionAttr(abc.ABC):
    """Descriptor base class"""
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
    def _validate(self):
        pass


class ScalarPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a scalar
    `~regions.PixCoord` object.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and value.isscalar):
            raise ValueError(f'The {self.name} must be a 0D PixCoord object')


class OneDPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`, which takes a one
    dimensional `regions.PixCoord` object.
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
            raise ValueError(f'The {self.name} must be a 0D SkyCoord object')


class OneDSky(RegionAttr):
    """
    Descriptor class for `~regions.SkyRegion`, which takes a one
    dimensional `~astropy.coordinates.SkyCoord` object.
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
            raise ValueError(f'The {self.name} must be a scalar astropy Quantity object')


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
    Descriptor class is for `~regions.CompoundSkyRegion`, which takes a
    `~regions.SkyRegion` object.
    """

    def _validate(self, value):
        if not isinstance(value, SkyRegion):
            raise ValueError(f'The {self.name} must be a SkyRegion object')


class Meta(dict):
    """
    A base class for region metadata.
    """

    def __init__(self, seq=None, **kwargs):
        super().__init__()

        if seq:
            if isinstance(seq, dict):
                for key, val in seq.items():
                    self.__setitem__(key, val)
            else:
                for key, val in seq:
                    self.__setitem__(key, val)

        if len(kwargs) > 0:
            for key, val in kwargs.items():
                self.__setitem__(key, val)

    def __setitem__(self, key, value):
        key = self.key_mapping.get(key, key)
        if key in self.valid_keys:
            super().__setitem__(key, value)
        else:
            raise KeyError(
                f"{key} is not a valid key for this class.")

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super().__getitem__(item)


class RegionMeta(Meta):
    """
    A dictionary subclass that holds the meta attributes of the region.
    """

    valid_keys = ['label', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'tag', 'comment', 'line', 'name', 'select',
                  'highlite', 'fixed', 'edit', 'move', 'rotate', 'delete',
                  'source', 'background', 'corr', 'type', 'text']

    key_mapping = {}


class RegionVisual(Meta):
    """
    A dictionary subclass which holds the visual attributes of the
    region.
    """

    valid_keys = ['color', 'dash', 'font', 'dashlist', 'symsize', 'symthick',
                  'symbol', 'fontsize', 'fontstyle', 'usetex', 'labelpos',
                  'labeloff', 'linewidth', 'linestyle', 'fill', 'line',
                  'textangle', 'fontweight']

    key_mapping = {'width': 'linewidth', 'point': 'symbol'}
