"""
Deal with numerical issues
"""

import numpy as np

def really_less(a, b, **kwargs):
    """
    Whether a < b, but with a small tolerance, i.e. whether a < b - eps for some eps > 0
    """
    
    # check whether a < b and a is not approximately b
    return (a < b) and not np.isclose(a, b, **kwargs)