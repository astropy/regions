# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the regions ASDF extension.
"""

import importlib
import re
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml
from asdf import get_config
from asdf.extension import Converter, ExtensionProxy

import regions.io.asdf.extension as extension_mod
from regions.io.asdf import converters
from regions.io.asdf.extension import (REGIONS_CONVERTERS,
                                       REGIONS_MANIFEST_URIS, get_extensions,
                                       get_resource_mappings)


def _manifest_tags():
    """
    Return the set of tag URIs defined by the regions manifests.
    """
    resource_manager = get_config().resource_manager
    tags = set()
    for manifest_uri in REGIONS_MANIFEST_URIS:
        manifest = yaml.safe_load(resource_manager[manifest_uri])
        tags.update(tag['tag_uri'] for tag in manifest['tags'])
    return tags


def _converter_tag_prefixes(converters_):
    """
    Return the set of unversioned tag URIs handled by ``converters_``.
    """
    return {tag.removesuffix('-*')
            for converter in converters_ for tag in converter.tags}


def test_converters_implement_asdf_interface():
    """
    Test that every registered converter implements the asdf ``Converter``
    interface.

    asdf rejects an entire extension that contains a converter which
    does not, so a single invalid converter would disable all regions
    ASDF support.
    """
    assert REGIONS_CONVERTERS
    for converter in REGIONS_CONVERTERS:
        assert isinstance(converter, Converter)


def test_extension_loads():
    """
    Test that the regions extension is accepted by asdf.
    """
    extensions = get_extensions()
    assert len(extensions) == len(REGIONS_MANIFEST_URIS)
    for extension in extensions:
        proxy = ExtensionProxy(extension)
        assert len(proxy.converters) == len(REGIONS_CONVERTERS)


def test_converter_tags_are_defined_by_the_manifest():
    """
    Test that every tag handled by a converter is defined by a manifest.
    """
    assert _converter_tag_prefixes(REGIONS_CONVERTERS)
    manifest_tags = _manifest_tags()
    for prefix in _converter_tag_prefixes(REGIONS_CONVERTERS):
        assert any(tag.startswith(prefix) for tag in manifest_tags), prefix


def test_manifest_tags_all_have_converters():
    """
    Test that every tag defined by a manifest has a converter.
    """
    prefixes = _converter_tag_prefixes(REGIONS_CONVERTERS)
    for tag in _manifest_tags():
        assert any(tag.startswith(prefix) for prefix in prefixes), tag


def test_resource_mappings_register_all_files():
    """
    Test that every schema and manifest file is registered with asdf and
    resolves to a resource with a matching ``id``.
    """
    resource_manager = get_config().resource_manager

    uris = set()
    for mapping in get_resource_mappings():
        uris.update(mapping)
    assert uris

    for uri in uris:
        assert uri in resource_manager
        assert yaml.safe_load(resource_manager[uri])['id'] == uri


def test_public_names_are_importable():
    """
    Test that every name in ``__all__`` is defined.

    The converters are imported conditionally, so ``__all__`` must be
    built to match.
    """
    assert converters.__all__
    for name in converters.__all__:
        assert hasattr(converters, name)


def test_documented_tags_match_the_manifest():
    """
    Test that the tags listed in the ASDF documentation match those
    defined by the manifests.

    The documentation table is written by hand, so it would otherwise
    go stale when a tag is added, removed, or renamed.
    """
    docs = (Path(__file__).parents[3] / 'docs' / 'user_guide' / 'asdf.rst')
    if not docs.is_file():
        pytest.skip('the documentation is not part of the installed package')

    documented = set(re.findall(r'``(regions/\w+/[\w.-]+)``',
                                docs.read_text()))
    assert documented == {tag.split(':')[-1] for tag in _manifest_tags()}


@pytest.mark.parametrize('missing', ['asdf'])
def test_missing_optional_dependency(missing):
    """
    Test the converter registration when an optional dependency is not
    installed.

    The modules are reloaded with ``find_spec`` patched to report the
    package as missing, and then reloaded again to restore them.
    """
    real_find_spec = importlib.util.find_spec

    def find_spec(name, *args, **kwargs):
        if name == missing:
            return None
        return real_find_spec(name, *args, **kwargs)

    try:
        with patch('importlib.util.find_spec', new=find_spec):
            importlib.reload(converters)
            if missing == 'asdf':
                assert not converters._ASDF_INSTALLED
                assert converters.__all__ == []
            else:
                assert converters._ASDF_INSTALLED
                assert converters.__all__
                importlib.reload(extension_mod)
                assert extension_mod.REGIONS_CONVERTERS
    finally:
        # Reload again to restore the original state
        importlib.reload(converters)
        importlib.reload(extension_mod)
