# Licensed under a 3-clause BSD style license - see LICENSE.rst

import asdf
import pytest
from numpy.testing import assert_array_equal


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
