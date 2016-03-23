import math
import numpy as np
from ..read_ds9 import parse_ds9

def test_physical():
    """
    should return a list of  region objects
    """
    r=parse_ds9("regions/ds9.physical.reg")
    assert (type(r) == type([]) )


def test_physical_span():
    """
    should return a list  region objects
    """
    r=parse_ds9("regions/ds9.physical.strip.reg")
    assert (type(r) == type([]) )

def test_physical_span():
    """
    should return a list  region objects
    """
    r=parse_ds9("regions/ds9.color.reg")
    assert (type(r) == type([]) )
