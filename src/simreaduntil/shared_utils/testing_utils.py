"""
Utilities to support testing
"""

from typing import Dict, Any
import numpy as np

def assert_dict_with_np_arrays_equal(d1: Dict[Any, np.ndarray], d2: Dict[Any, np.ndarray]):
    """
    Assert that two dicts of np arrays are equal
    
    Args:
        d1: dict of np arrays
        d2: dict of np arrays
    """
    assert d1.keys() == d2.keys()
    for k in d1.keys():
        assert np.array_equal(d1[k], d2[k])