# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module proves classes to handle region metadata.
"""

__all__ = ['Meta', 'RegionMeta', 'RegionVisual']


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
