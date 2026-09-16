# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
ASDF converter for PixCoord object.
"""
from asdf.extension import Converter

__all__ = [
    'PixCoordConverter',
]


class PixCoordConverter(Converter):
    """
    ASDF converter for pixel coordinates.
    """

    tags = ('tag:astropy.org:regions/pixcoord-*',)
    types = ('regions.core.pixcoord.PixCoord',)

    def to_yaml_tree(self, obj, tag, ctx):  # noqa: ARG002

        return {
            'x': obj.x,
            'y': obj.y,
        }

    def from_yaml_tree(self, node, tag, ctx):  # noqa: ARG002
        from regions.core.pixcoord import PixCoord
        return PixCoord(
            x=node['x'],
            y=node['y'],
        )
