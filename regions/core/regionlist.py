# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides a RegionList class.
"""

from .core import Region
from .registry import RegionsRegistry

__all__ = ['RegionList']


class RegionList:
    """
    Class to hold a list of `~regions.Region` objects.

    Parameters
    ----------
    regions : list of `~region.Region`
        The list of region objects.
    """

    def __init__(self, regions):
        self.regions = regions

    def __getitem__(self, index):
        newregions = self.regions[index]
        if isinstance(newregions, Region):  # one item
            return newregions
        else:
            newcls = object.__new__(self.__class__)
            newcls.regions = newregions
            return newcls

    def __len__(self):
        return len(self.regions)

    def append(self, item):
        self.regions.append(item)

    def extend(self, item):
        self.regions.extend(item)

    def insert(self, index, item):
        self.regions.insert(index, item)

    def reverse(self):
        self.regions.reverse()

    def pop(self, index=-1):
        return self.regions.pop(index)

    def copy(self):
        newcls = object.__new__(self.__class__)
        newcls.regions = self.regions.copy()
        return newcls

    @classmethod
    def read(cls, filename, format=None, cache=False, **kwargs):
        """
        Read in a regions file.
        """
        return RegionsRegistry.read(filename, cls.__name__, format=format,
                                    cache=cache, **kwargs)

    @classmethod
    def parse(cls, data, format=None, **kwargs):
        """
        Parse a regions string or table.
        """
        return RegionsRegistry.parse(data, cls.__name__, format=format,
                                     **kwargs)

    def write(self, filename, format=None, **kwargs):
        """
        Write to a regions file.
        """
        return RegionsRegistry.write(self.regions, filename,
                                     self.__class__.__name__, format=format,
                                     **kwargs)


# reuse the `list` API docs for these methods
methods = ('append', 'extend', 'insert', 'reverse', 'pop', 'copy')
for method in methods:
    getattr(RegionList, method).__doc__ = getattr(list, method).__doc__
