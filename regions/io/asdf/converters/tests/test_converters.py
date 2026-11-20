# Licensed under a 3-clause BSD style license - see LICENSE.rst

import asdf
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from . import examples


@pytest.mark.parametrize('example', [examples.pixcoord,
                                     examples.pixcoord_array],
                         ids=lambda example: example.__name__)
def test_region_converters(tmp_path, example):
    """
    Test that the Region converters can round-trip a Region object.
    """
    region, pars = example()
    filename = tmp_path / 'region.asdf'
    with asdf.AsdfFile() as af:
        af['region'] = region
        af.write_to(filename)

    with asdf.open(filename) as af:
        region2 = af['region']

    assert isinstance(region2, type(region))

    for name in pars:
        value, value2 = getattr(region, name), getattr(region2, name)
        assert np.shape(value2) == np.shape(value)
        assert np.asarray(value2).dtype == np.asarray(value).dtype
        assert_array_equal(value2, value)
