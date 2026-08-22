# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from regions.io.asdf.converters.tests import examples


@pytest.fixture
def test_pixcoord():
    """
    Test that PixCoord can be serialized and deserialized correctly.
    """
    return examples.pixcoord()
