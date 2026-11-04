# Licensed under a 3-clause BSD style license - see LICENSE.rst

import asdf
import pytest
from numpy.testing import assert_array_equal, assert_equal


@pytest.fixture
def obj(request):
    """
    A pytest fixture that returns a 'Region' instance and the
    list of parameters to test.
    """
    return request.getfixturevalue(request.param)


region_params = pytest.mark.parametrize('obj', [
    'test_pixcoord',
    'test_pixcoord_array',
], indirect=True)


region_meta = pytest.mark.parametrize('obj', [
    'test_meta',
    'test_visual_meta',
], indirect=True)


@region_params
def test_region_converters(tmp_path, obj):
    """
    Test that the Region converters can round-trip a Region object.
    """
    region, pars = obj

    with asdf.AsdfFile() as af:
        af['region'] = region
        af.write_to(tmp_path / 'region.asdf')

        with asdf.open(tmp_path / 'region.asdf') as af:
            region2 = af['region']
            for parameter in pars:
                assert_array_equal(getattr(region, parameter),
                                   getattr(region2, parameter))
    if 'meta' in pars:
        # exclude PixCoord from this test.
        assert_equal(region.meta, region2.meta)
        assert_equal(region.visual, region2.visual)


@region_meta
def test_meta(tmp_path, obj):
    """
    Test metadata can be serialized and deserialized correctly.
    """
    metadata, params = obj
    with asdf.AsdfFile() as af:
        af['meta'] = metadata
        af.write_to(tmp_path / 'meta.asdf')

        with asdf.open(tmp_path / 'meta.asdf') as af:
            metadata2 = af['meta']
            assert_equal(metadata, metadata2)
