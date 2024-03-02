import numpy as np
import pytest
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked

from simreaduntil.simulator.simulator_params import SimParams

def test_sim_params():
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=1.2, long_gap_length=1.2, prob_long_gap=0.15, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=10, min_chunk_size=4, default_unblock_duration=1.2, seed=0,
    )
    
    # set to some random values
    sim_params.set(bp_per_second=1000, default_unblock_duration=0.2, min_chunk_size=100, seed=2)
    assert sim_params.bp_per_second == 1000
    assert sim_params.default_unblock_duration == 0.2
    assert sim_params.min_chunk_size == 100
    assert sim_params._initial_seed == 2
    assert sim_params.gap_samplers["channel_0"].short_gap_length == 1.2
    
    print(sim_params)