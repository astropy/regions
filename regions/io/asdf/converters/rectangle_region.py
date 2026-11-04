# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
ASDF converter for Rectangle regions in pixel and sky coordinates.
"""
from asdf.extension import Converter

__all__ = [
    'RectangleConverter',
]


class RectangleConverter(Converter):
    """
    ASDF converter for rectangle regions in pixel and sky coordinates.
    """

    tags = ('tag:astropy.org:regions/rectangle_region-*',)
    types = ('regions.shapes.rectangle.RectanglePixelRegion',
             'regions.shapes.rectangle.RectangleSkyRegion')

    def to_yaml_tree(self, obj, tag, ctx):  # noqa: ARG002

        node = {
            'center': obj.center,
            'width': obj.width,
            'height': obj.height,
            'meta': obj.meta,
            'visual': obj.visual}
        return node

    def from_yaml_tree(self, node, tag, ctx):  # noqa: ARG002
        from regions.core.pixcoord import PixCoord
        center = node['center']
        width = node['width']
        height = node['height']
        meta = node['meta']
        visual = node['visual']
        if isinstance(center, PixCoord):
            from regions.shapes.rectangle import RectanglePixelRegion
            return RectanglePixelRegion(center=center,
                                         width=width,
                                         height=height,
                                         meta=meta, visual=visual)
        else:
            from regions.shapes.rectangle import RectangleSkyRegion
            return RectangleSkyRegion(center=center,
                                       width=width,
                                       height=height,
                                       meta=meta, visual=visual)
