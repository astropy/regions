# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
ASDF extension for regions.
"""
import importlib.resources as importlib_resources

from asdf.extension import ManifestExtension
from asdf.resource import DirectoryResourceMapping

from regions.io.asdf.converters import circle_region
from regions.io.asdf.converters import metadata
from regions.io.asdf.converters import pixel_coords
from regions.shapes.circle import CircleSphericalSkyRegion


__all__ = ['REGIONS_CONVERTERS']


REGIONS_CONVERTERS = [
    pixel_coords.PixCoordConverter(),
    circle_region.CircleConverter(),
    # circle.CircleAnnulusPixelRegion(),
    # circle.CircleAnnulusSkyRegion(),
    # circle.CircleAnnulusSphericalSkyRegion(),
    metadata.RegionMetaConverter(),
    metadata.RegionVisualConverter(),
]

# The order here is important; asdf will prefer to use extensions
# that occur earlier in the list.
REGIONS_MANIFEST_URIS = [
    'asdf://astropy.org/regions/manifests/regions-1.0.0',
]


def get_extensions():
    """
    Get the regions extension.
    This method is registered with the asdf.extensions entry point.

    Returns
    -------
    list
        A list of ASDF extensions.
    """
    return [
        ManifestExtension.from_uri(
            uri,
            converters=REGIONS_CONVERTERS,
        )
        for uri in REGIONS_MANIFEST_URIS
    ]


def get_resource_mappings():
    """
    Get the resource mapping instances for the regions schemas
    and manifests.  This method is registered with the
    asdf.resource_mappings entry point.

    Returns
    -------
    list
        A list of collections.abc.Mapping of ASDF resource mappings.
    """
    from regions.io.asdf import resources

    resources_root = importlib_resources.files(resources)

    return [
        DirectoryResourceMapping(resources_root / 'schemas',
                                 'asdf://astropy.org/regions/schemas',
                                 recursive=True),
        DirectoryResourceMapping(resources_root / 'manifests',
                                 'asdf://astropy.org/regions/manifests'),
    ]
