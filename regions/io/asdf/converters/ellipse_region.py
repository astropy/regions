# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
ASDF converter for Ellipse regions in pixel and sky coordinates.
"""
from asdf.extension import Converter

__all__ = [
    'EllipseConverter',
]


class EllipseConverter(Converter):
    """
    ASDF converter for ellipse regions in pixel and sky coordinates.
    """

    tags = ('tag:astropy.org:regions/ellipse_region-*',)
    types = ('regions.shapes.ellipse.EllipsePixelRegion',
             'regions.shapes.ellipse.EllipseSkyRegion')

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
            from regions.shapes.ellipse import EllipsePixelRegion
            return EllipsePixelRegion(center=center,
                                      width=width,
                                      height=height,
                                      meta=meta, visual=visual)
        else:
            from regions.shapes.ellipse import EllipseSkyRegion
            return EllipseSkyRegion(center=center,
                                    width=width,
                                    height=height,
                                    meta=meta, visual=visual)
