# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides the astropy unified I/O layer to a RegionList
class.
"""

from astropy.io.registry import UnifiedReadWrite, read, write


__all__ = []


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
