import math
import numpy as np
from ..read_ds9 import parse_ds9

def test_physical():
    """
    should return a list circle region objects
    """
    r=parse_ds9("regions/ds9.physical.reg")
    assert (type(r) == type([]) )
    
