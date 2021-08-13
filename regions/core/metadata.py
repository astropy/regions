# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module proves classes to handle region metadata.
"""

__all__ = ['Meta', 'RegionMeta', 'RegionVisual']


class Meta(dict):
    """
    A base class for region metadata.
    """

    valid_keys = []
    key_mapping = {}

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
            raise KeyError(f'{key} is not a valid key for this class.')

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super().__getitem__(item)


class RegionMeta(Meta):
    """
    A dictionary subclass that holds the meta attributes of the region.
    """

    valid_keys = ['background', 'comment', 'corr', 'delete', 'edit', 'fixed',
                  'frame', 'highlite', 'include', 'label', 'line', 'move',
                  'name', 'range', 'restfreq', 'rotate', 'select', 'source',
                  'tag', 'text', 'type', 'veltype']

    key_mapping = {}


class RegionVisual(Meta):
    """
    A dictionary subclass which holds the visual attributes of the
    region.
    """

    valid_keys = ['color', 'dash', 'dashlist', 'fill', 'font', 'fontsize',
                  'fontstyle', 'fontweight', 'labeloff', 'labelpos', 'line',
                  'linestyle', 'linewidth', 'symbol', 'symsize', 'symthick',
                  'textangle', 'usetex', 'default_style']

    key_mapping = {'point': 'symbol', 'width': 'linewidth'}

    def _define_default_mpl_kwargs(self, artist):
        """
        Define the default matplotlib kwargs for the specified artist.

        The kwargs depend on the value of self.visual['default_style'],
        which can be set when reading region files. If this is not set
        or set to `mpl' or `None`, then the matplotlib defaults will be
        used, with the exception that Patch and Line2D fill is turned
        off.

        Parameters
        ----------
        artist : {'Text', 'Line2D', 'Patch'}
            The matplotlib artist type.

        Returns
        -------
        result : dict
            A dictionary of matplotlib keyword arguments.
        """
        kwargs = {}

        default_style = self.get('default_style', None)
        if default_style is None or default_style == 'mpl':
            # do not fill by default, which is the only change from
            # matplotlib defaults
            if artist == 'Patch':
                kwargs['fill'] = False
            elif artist == 'Line2D':
                kwargs['fillstyle'] = 'none'
                kwargs['marker'] = 'o'
            return kwargs

        # 'ds9' style is set when reading from ds9 region files
        elif default_style == 'ds9':
            kwargs['color'] = '#00ff00'  # green
            if artist == 'Text':
                kwargs['ha'] = 'center'  # text horizontal alignment
                kwargs['va'] = 'center'  # text vertical alignment
            elif artist == 'Line2D':
                from ..io.ds9.core import valid_symbols_ds9

                kwargs['marker'] = valid_symbols_ds9['boxcircle']
                kwargs['markersize'] = 11
                kwargs['markeredgecolor'] = kwargs['color']
                kwargs['fillstyle'] = 'none'
            elif artist == 'Patch':
                kwargs['edgecolor'] = kwargs.pop('color')
                kwargs['fill'] = False

            else:
                raise ValueError('invalid visual["default"] value')

        return kwargs
