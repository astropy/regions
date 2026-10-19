"""
ASDF converter for Circle regions in pixel coordinates.
"""
from asdf.extension import Converter
from regions.shapes.circle import CircleSphericalSkyRegion

__all__ = [
    'CircleConverter',
]


class CircleConverter(Converter):
    """
    ASDF converter for circle regions in pixel coordinates.
    """

    tags = ('tag:astropy.org:regions/circle_region-*',)
    types = ('regions.shapes.circle.CirclePixelRegion',
             'regions.shapes.circle.CircleSkyRegion',
             'regions.shapes.circle.CircleSphericalSkyRegion',)


    def to_yaml_tree(self, obj, tag, ctx):  # noqa: ARG002
        from regions.shapes.circle import CircleSphericalSkyRegion
        spherical = False
        if isinstance (obj, CircleSphericalSkyRegion):
            spherical = True

        return {
            'center': obj.center,
            'radius': obj.radius,
            'spherical': spherical
        }

    def from_yaml_tree(self, node, tag, ctx):  # noqa: ARG002
        from regions.core.pixcoord import PixCoord

        center=node['center'],
        radius=node['radius'],
        spherical=node['spherical']
        if isinstance(center, PixCoord):
            from regions.shapes.circle import CirclePixelRegion
            return CirclePixelRegion(center=center, radius=radius)
        if spherical:
            from regions.shapes.circle import CircleSphericalSkyRegion
            return CircleSphericalSkyRegion(center=center, radius=radius)
        else:
            from regions.shapes.circle import CircleSkyRegion
            return CircleSkyRegion(center=center, radius=radius)
