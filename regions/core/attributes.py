# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
import six
import weakref

from astropy.coordinates import SkyCoord
from astropy.units import Quantity

import numpy as np

from .pixcoord import PixCoord
from .core import PixelRegion, SkyRegion

__all__ = ['RegionMeta', 'RegionVisual']

"""
Several custom descriptor classes are present here for attribute validation of
region classes.
Also, contains RegionMeta and RegionVisual classes to handle meta data of regions.
"""


@six.add_metaclass(abc.ABCMeta)
class RegionAttr(object):

    def __init__(self, name):
        self._name = name
        self._values = weakref.WeakKeyDictionary()

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return self._values.get(instance, None)

    def __set__(self, instance, value):
        self._validate(value)
        self._values[instance] = value

    def _validate(self, value):
        raise NotImplementedError


class ScalarPix(RegionAttr):

    def _validate(self, value):
        if not(isinstance(value, PixCoord) and value.isscalar):
            raise ValueError('The {} must be a 0D PixCoord object'
                             .format(self._name))


class OneDPix(RegionAttr):

    def _validate(self, value):
        if not(isinstance(value, PixCoord) and not value.isscalar
               and value.x.ndim == 1):
            raise ValueError('The {} must be a 1D PixCoord object'
                             .format(self._name))


class AnnulusCenterPix(object):

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg1 = getattr(instance, 'region1')
        return getattr(reg1, 'center')

    def __set__(self, instance, value):

        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if isinstance(value, PixCoord) and value.isscalar:
            setattr(reg1, 'center', value)
            setattr(reg2, 'center', value)
        else:
            raise ValueError('The center must be a 0D PixCoord object')


class ScalarLength(RegionAttr):

    def _validate(self, value):
        if not np.isscalar(value):
            raise ValueError(
                'The {} must be a scalar numpy/python number'.format(self._name))


class AnnulusInnerScalarLength(object):

    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg1 = getattr(instance, 'region1')
        return getattr(reg1, self._name)

    def __set__(self, instance, value):
        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if np.isscalar(value):
            if getattr(reg2, self._name) < value:
                raise ValueError("The inner {0} must be less than the outer {0}"
                                .format(self._name)
                             )
            else:
                setattr(reg1, self._name, value)
        else:
            raise ValueError('The inner {} must be a scalar numpy/python number'
                             .format(self._name))


class AnnulusOuterScalarLength(object):

    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg2 = getattr(instance, 'region2')
        return getattr(reg2, self._name)

    def __set__(self, instance, value):
        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if np.isscalar(value):
            if getattr(reg1, self._name) > value:
                raise ValueError("The outer {0} must be greater than the outer"
                                 " {0}".format(self._name)
                                 )
            else:
                setattr(reg2, self._name, value)
        else:
            raise ValueError('The outer {} must be a scalar numpy/python number'
                             .format(self._name))


class ScalarSky(RegionAttr):

    def _validate(self, value):
        if not(isinstance(value, SkyCoord) and value.isscalar):
            raise ValueError('The {} must be a 0D SkyCoord object'.
                             format(self._name))


class OneDSky(RegionAttr):

    def _validate(self, value):
        if not(isinstance(value, SkyCoord) and value.ndim == 1):
            raise ValueError('The {} must be a 1D SkyCoord object'.
                             format(self._name))


class AnnulusCenterSky(object):

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg1 = getattr(instance, 'region1')
        return getattr(reg1, 'center')

    def __set__(self, instance, value):

        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if isinstance(value, SkyCoord) and value.isscalar:
            setattr(reg1, 'center', value)
            setattr(reg2, 'center', value)
        else:
            raise ValueError('The center must be a 0D SkyCoord object')


class QuantityLength(RegionAttr):

    def _validate(self, value):
        if not(isinstance(value, Quantity) and value.isscalar):
            raise ValueError('The {} must be a scalar astropy Quantity object'
                             .format(self._name))


class AnnulusInnerQuantityLength(object):

    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg1 = getattr(instance, 'region1')
        return getattr(reg1, self._name)

    def __set__(self, instance, value):
        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if isinstance(value, Quantity) and value.isscalar:
            if getattr(reg2, self._name) < value:
                raise ValueError("The inner {0} must be less than the outer {0}"
                                 .format(self._name)
                                 )
            else:
                setattr(reg1, self._name, value)
        else:
            raise ValueError('The inner {} must be a scalar astropy Quantity '
                             'object'.format(self._name))


class AnnulusOuterQuantityLength(object):

    def __init__(self, name):
        self._name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg2 = getattr(instance, 'region2')
        return getattr(reg2, self._name)

    def __set__(self, instance, value):
        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if isinstance(value, Quantity) and value.isscalar:
            if getattr(reg1, self._name) > value:
                raise ValueError("The inner {0} must be less than the outer {0}"
                                 .format(self._name)
                                 )
            else:
                setattr(reg2, self._name, value)
        else:
            raise ValueError('The outer {} must be a scalar astropy Quantity '
                             'object'.format(self._name))


class AnnulusAngle(object):

    def __get__(self, instance, owner):
        if instance is None:
            return self
        reg1 = getattr(instance, 'region1')
        return getattr(reg1, 'angle')

    def __set__(self, instance, value):

        reg1 = getattr(instance, 'region1')
        reg2 = getattr(instance, 'region2')

        if isinstance(value, Quantity) and value.isscalar:
            setattr(reg1, 'angle', value)
            setattr(reg2, 'angle', value)
        else:
            raise ValueError('The angle must be a scalar astropy quantity object')


class CompoundRegionPix(RegionAttr):

    def _validate(self, value):
        if not isinstance(value, PixelRegion):
            raise ValueError('The {} must be a PixelRegion object'
                             .format(self._name))


class CompoundRegionSky(RegionAttr):

    def _validate(self, value):
        if not isinstance(value, SkyRegion):
            raise ValueError('The {} must be a SkyRegion object'
                             .format(self._name))


class RegionMeta(dict):
    """
    A python dictionary subclass which holds the meta attributes of the region.
    """
    valid_keys = ['label', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'tag', 'comment', 'line', 'name',
                  'select', 'highlite', 'fixed', 'edit', 'move', 'rotate',
                  'delete', 'source', 'background', 'corr', 'type'
                  ]

    key_mapping = {}

    def __setitem__(self, key, value):
        key = self.key_mapping.get(key, key)
        if key in self.valid_keys:
            super(RegionMeta, self).__setitem__(key, value)
        else:
            raise KeyError("{} is not a valid meta key for region.".format(key))

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super(RegionMeta, self).__getitem__(item)


class RegionVisual(dict):
    """
    A python dictionary subclass which holds the visual attributes of the region.
    """
    valid_keys = ['color', 'dash', 'font', 'dashlist', 'symsize', 'symthick',
                  'symbol', 'fontsize', 'fontstyle', 'usetex', 'labelpos',
                  'labeloff', 'linewidth', 'linestyle', 'fill', 'line',
                  'textangle', 'fontweight']

    key_mapping = {'width': 'linewidth', 'point': 'symbol'}

    def __setitem__(self, key, value):
        key = self.key_mapping.get(key, key)
        if key in self.valid_keys:
            super(RegionVisual, self).__setitem__(key, value)
        else:
            raise KeyError("{} is not a valid visual meta key for region.".format(key))

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super(RegionVisual, self).__getitem__(item)