# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides a RegionList class.
"""

from astropy.io.registry import (UnifiedReadWrite, UnifiedReadWriteMethod,
                                 read, write)

from .core import Region

__all__ = ['RegionList']


class RegionListRead(UnifiedReadWrite):
    """
    Read and parse a region file and return as a `~regions.RegionList`.

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data reader. If supplied
        the first argument is typically the input filename.

    **kwargs : dict, optional
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~regions.RegionList`
        A `~regions.RegionList` corresponding to the region file
        contents.
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'read')

    def __call__(self, *args, **kwargs):
        return read(self._cls, *args, **kwargs)


class RegionListWrite(UnifiedReadWrite):
    """
    Write the `~regions.RegionList` object to a file in the specified
    format.

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied,
        the first argument is the output filename.

    **kwargs : dict, optional
        Keyword arguments passed through to data writer.
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'write')

    def __call__(self, *args, **kwargs):
        write(self._instance, *args, **kwargs)


class RegionList:
    """
    Class to hold a list of `~regions.Region` objects.

    Parameters
    ----------
    regions : list of `~region.Region`
        The list of region objects.
    """

    # Unified I/O read and write methods
    read = UnifiedReadWriteMethod(RegionListRead)
    write = UnifiedReadWriteMethod(RegionListWrite)

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


# use list API docs for these methods
methods = ('append', 'extend', 'insert', 'reverse', 'pop', 'copy')
for method in methods:
    getattr(RegionList, method).__doc__ = getattr(list, method).__doc__
