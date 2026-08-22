# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
ASDF extension for regions.
"""
import importlib.resources as importlib_resources

from asdf.extension import ManifestExtension
from asdf.resource import DirectoryResourceMapping

from regions.io.asdf.converters.pixel_coords import PixCoordConverter

__all__ = [
    'REGIONS_CONVERTERS',
    'REGIONS_MANIFEST_URIS',
    'get_extensions',
    'get_resource_mappings',
]


REGIONS_CONVERTERS = [
    PixCoordConverter(),
]

# The order here is important. asdf will prefer to use extensions
# that occur earlier in the list.
REGIONS_MANIFEST_URIS = [
    'asdf://astropy.org/regions/manifests/regions-1.0.0',
]


def get_extensions():
    """
    Get the regions ASDF extension.

    This function is registered with the ``asdf.extensions`` entry point.

    Returns
    -------
    extensions : list of ``asdf.extension.ManifestExtension``
        The ASDF extensions, one for each regions manifest.
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
    and manifests.

    This function is registered with the ``asdf.resource_mappings``
    entry point.

    Returns
    -------
    directory_resource_mapping : list of asdf.resource.DirectoryResourceMapping
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
