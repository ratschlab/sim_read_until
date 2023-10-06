
import numpy as np
import pandas as pd
from pytest import approx
from simreaduntil.simulator.gap_sampling.gap_sampler_per_window_until_blocked import compute_time_and_aggregation_windows

def test_compute_time_and_aggregation_windows():
    assert compute_time_and_aggregation_windows(100, 4, 0.5) == (
        approx(np.array([[0, 25], [25, 50], [50, 75], [75, 100]])),
        approx(np.array([[0, 37.5], [12.5, 62.5], [37.5, 87.5], [62.5, 100]]))
    )

    assert compute_time_and_aggregation_windows(100, 1, 0.5) == (approx(np.array([[0, 100]])), approx(np.array([[0, 100]])))
    assert compute_time_and_aggregation_windows(100, 1, 1) == (
        approx(np.array([[0, 100]])),
        approx(np.array([[0, 100]]))
    )

    assert compute_time_and_aggregation_windows(100, 2, 0.5) == (
        approx(np.array([[0, 50], [50, 100]])),
        approx(np.array([[0, 75], [25, 100]]))
    )
    assert compute_time_and_aggregation_windows(100, 2, 0) == (
        approx(np.array([[0, 50], [50, 100]])),
        approx(np.array([[0, 50], [50, 100]]))
    )
    assert compute_time_and_aggregation_windows(100, 2, 1) == (
        approx(np.array([[0, 50], [50, 100]])),
        approx(np.array([[0, 100], [0, 100]]))
    )
    