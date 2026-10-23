"""
ASDF converter for RegionMeta and RegionVisual.
"""
from asdf.extension import Converter
from regions.io.asdf.converters.tests.test import obj


__all__ = [
    'RegionMetaConverter',
    'RegionVisualConverter',
]


class RegionMetaConverter(Converter):
    """
    ASDF converter for RegionMeta.
    """

    tags = ('tag:astropy.org:regions/meta-*',)
    types = ('regions.core.metadata.RegionMeta',)


    def to_yaml_tree(self, obj, tag, ctx):  # noqa: ARG002

        return {'meta': dict(obj.items())}

    def from_yaml_tree(self, node, tag, ctx):  # noqa: ARG002
        from regions.core.metadata import RegionMeta
        return RegionMeta(**node['meta'])


class RegionVisualConverter(Converter):
    """
    ASDF converter for RegionVisual.
    """

    tags = ('tag:astropy.org:regions/visual-*',)
    types = ('regions.core.metadata.RegionVisual',)


    def to_yaml_tree(self, obj, tag, ctx):  # noqa: ARG002
        return {'visual_meta': dict(obj.items())}

    def from_yaml_tree(self, node, tag, ctx):  # noqa: ARG002
        from regions.core.metadata import RegionVisual
        return RegionVisual(**node['visual_meta'])
    