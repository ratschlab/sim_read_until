import numpy as np
import pytest

from simreaduntil.shared_utils.testing_utils import assert_dict_with_np_arrays_equal

def test_check_dict_with_np_arrays_equal():
    assert_dict_with_np_arrays_equal({"a": np.array([1, 2, 3])}, {"a": np.array([1, 2, 3])})
    
    with pytest.raises(AssertionError):
        assert_dict_with_np_arrays_equal({"a": np.array([1, 2, 3])}, {"a": np.array([1, 2, 3]), "b": np.array([4, 5])})
        
    with pytest.raises(AssertionError):
        assert_dict_with_np_arrays_equal({"a": np.array([1])}, {"a": np.array([1, 2, 3])})