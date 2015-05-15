import math
import numpy as np
from ..circle import CircularRegion

def test_basic():
    c = CircularRegion((3,4), 2)
    np.testing.assert_allclose(c.area, 4 * math.pi)
    assert (3.4, 4.1) in c
    assert not (10., 3.) in c