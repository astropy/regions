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
    """
     Meta descriptor class for attribute of an `~regions.Region`
     which makes sure that it's value is valid all the time.
    """

    def __init__(self, name):
        self._name = name

        # WeakKeyDictionary has an object instance as the key
        # and the key value pair remains until the object instance is not
        # free from memory, thus prevents memory leak.
        self._values = weakref.WeakKeyDictionary()

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return self._values.get(instance, None)

    def __set__(self, instance, value):
        self._validate(value)
        self._values[instance] = value

    def _validate(self, value):
        """
        Validates the value of the attribute. Raises an exception if invalid
        else does nothing.
        """
        raise NotImplementedError


class ScalarPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`  which allows values
    to be a scalar `~regions.PixCoord` object.
    """

    def _validate(self, value):
        if not(isinstance(value, PixCoord) and value.isscalar):
            raise ValueError('The {} must be a 0D PixCoord object'
                             .format(self._name))


class OneDPix(RegionAttr):
    """
    Descriptor class for `~regions.PixelRegion`  which allows values
    to be a one dimensional `regions.PixCoord` object.
    """

    def _validate(self, value):
        if not(isinstance(value, PixCoord) and not value.isscalar
               and value.x.ndim == 1):
            raise ValueError('The {} must be a 1D PixCoord object'
                             .format(self._name))


class AnnulusCenterPix(object):
    """
    This descriptor class is for the ``center`` of an
    ``annulus`` `~regions.PixelRegion`. It allows the value to be a scalar
    `~regions.PixCoord` object. It also makes sure that ``region1`` and
    ``region2`` are in sync in case of an updation.
    """

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
    """
    Descriptor class for `~regions.PixelRegion`  which allows
    values to be a scalar python/numpy number.
    """

    def _validate(self, value):
        if not np.isscalar(value):
            raise ValueError(
                'The {} must be a scalar numpy/python number'.format(self._name))


class AnnulusInnerScalarLength(object):
    """
    This descriptor class is for an inner length of an
    ``annulus`` `~regions.PixelRegion`. It allows the value to be a scalar
    python/numpy number and makes sure that it is less than the outer
    length of the annulus.
    """

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
    """
    This descriptor class is for an outer length of an
    ``annulus`` `~regions.PixelRegion`. It allows the values to be a scalar
    python/numpy number and makes sure that it is greater than the inner
    length of the annulus.
    """

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
    """
    Descriptor class for `~regions.SkyRegion` which allows values to be a
    scalar `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not(isinstance(value, SkyCoord) and value.isscalar):
            raise ValueError('The {} must be a 0D SkyCoord object'.
                             format(self._name))


class OneDSky(RegionAttr):
    """
    Descriptor class for `~regions.SkyRegion` which allows values to be a
    one dimensional `~astropy.coordinates.SkyCoord` object.
    """

    def _validate(self, value):
        if not(isinstance(value, SkyCoord) and value.ndim == 1):
            raise ValueError('The {} must be a 1D SkyCoord object'.
                             format(self._name))


class AnnulusCenterSky(object):
    """
    This descriptor class is for the ``center`` of an
    ``annulus`` `~regions.SkyRegion`. It allows the value to be a scalar
    `~astropy.coordinates.SkyCoord` object. It also makes sure that ``region1``
    and ``region2`` are in sync in case of an updation.
    """

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
    """
    Descriptor class for `~regions.SkyRegion`  which allows
    values to be a scalar `~astropy.units.Quantity` object.
    """

    def _validate(self, value):
        if not(isinstance(value, Quantity) and value.isscalar):
            raise ValueError('The {} must be a scalar astropy Quantity object'
                             .format(self._name))


class AnnulusInnerQuantityLength(object):
    """
    This descriptor class is for an inner length of an
    ``annulus`` `~regions.SkyRegion`. It allows the value to be a scalar
    `astropy.units.Quantity` object and makes sure that it is less than the outer
    length of the annulus.
    """

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
    """
    This descriptor class is for an outer length of an
    ``annulus`` `~regions.SkyRegion`. It allows the value to be a scalar
    `astropy.units.Quantity` object and makes sure that it is less than the
    inner length of the annulus.
    """

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
    """
    This descriptor class is for the ``center`` of an
    ``annulus`` `~regions.SkyRegion`. It allows the value to be a scalar
    `~astropy.units.Quantity` object. It also makes sure that ``region1``
    and ``region2`` are in sync in case of an updation.
    """

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
    """
    Descriptor class for `~regions.CompoundPixelRegion` which allows values
    to be `~regions.PixelRegion` object.
    """

    def _validate(self, value):
        if not isinstance(value, PixelRegion):
            raise ValueError('The {} must be a PixelRegion object'
                             .format(self._name))


class CompoundRegionSky(RegionAttr):
    """
    Descriptor class is for `~regions.CompoundSkyRegion` which allows values
    to be `~regions.SkyRegion` object.
    """

    def _validate(self, value):
        if not isinstance(value, SkyRegion):
            raise ValueError('The {} must be a SkyRegion object'
                             .format(self._name))


@six.add_metaclass(abc.ABCMeta)
class Meta(dict):

    def __init__(self, seq=None, **kwargs):

        super(Meta, self).__init__()

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
            super(Meta, self).__setitem__(key, value)
        else:
            raise KeyError(
                "{} is not a valid key for this class.".format(key))

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super(Meta, self).__getitem__(item)


class RegionMeta(Meta):
    """
    A python dictionary subclass which holds the meta attributes of the region.
    """
    valid_keys = ['label', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'tag', 'comment', 'line', 'name',
                  'select', 'highlite', 'fixed', 'edit', 'move', 'rotate',
                  'delete', 'source', 'background', 'corr', 'type'
                  ]

    key_mapping = {}


class RegionVisual(Meta):
    """
    A python dictionary subclass which holds the visual attributes of the region.
    """
    valid_keys = ['color', 'dash', 'font', 'dashlist', 'symsize', 'symthick',
                  'symbol', 'fontsize', 'fontstyle', 'usetex', 'labelpos',
                  'labeloff', 'linewidth', 'linestyle', 'fill', 'line',
                  'textangle', 'fontweight']

    key_mapping = {'width': 'linewidth', 'point': 'symbol'}
