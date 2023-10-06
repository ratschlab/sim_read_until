

from simreaduntil.shared_utils.numerics import really_less


def test_really_less():
    # default values for np.isclose: rtol=1e-5, atol=1e-8
    assert really_less(1 - 0.01, 1)
    assert not really_less(1 + 0.01, 1)
    assert not really_less(1, 1)
    
    assert not really_less(1 + 1e-16, 1)
    
    assert not really_less(1 - 1e-16, 1)
    assert not really_less(1 - 1e-5, 1)
    assert really_less(1 - 1e-4, 1, rtol=1e-5, atol=1e-8)