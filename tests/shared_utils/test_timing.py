import time
import numpy as np
import pytest
from pytest import approx

from simreaduntil.shared_utils.timing import ElapsedTimer, StepAndTotalElapsedTimer, cur_ns_time
from simreaduntil.shared_utils.utils import is_sorted
from simreaduntil.simulator.utils import in_interval

# Note: adding debugging breakpoints here distorts the times and thus some assertions may fail
    
def test_elapsed_timer():
    timer = ElapsedTimer()
    assert timer() < 0.01
    time.sleep(0.2)
    assert in_interval(timer(), (0.2, 0.22))
    time.sleep(0.3)
    assert in_interval(timer(reset=False), (0.3, 0.32))
    time.sleep(0.2)
    assert in_interval(timer(), (0.5, 0.52))
    
def test_step_and_total_elapsed_timer():
    timer = StepAndTotalElapsedTimer()
    
    assert timer(reset=True, do_print=False)[0] < 1e-4
    time.sleep(0.2)
    time_since_last_call, time_since_last_reset = timer(reset=True, do_print=False)
    assert time_since_last_call == approx(time_since_last_reset)
    assert in_interval(time_since_last_call, (0.2, 0.22))
    assert in_interval(time_since_last_reset, (0.2, 0.22))
    
    timer(reset=True, do_print=False)
    time.sleep(0.3)
    time_since_last_call, time_since_last_reset = timer(reset=False, do_print=False)
    assert in_interval(time_since_last_call, (0.3, 0.32))
    assert in_interval(time_since_last_reset, (0.3, 0.32))
    
    time.sleep(0.2)
    time_since_last_call, time_since_last_reset = timer(do_print=False)
    assert in_interval(time_since_last_call, (0.2, 0.22))
    assert in_interval(time_since_last_reset, (0.5, 0.52))
    
    timer(reset=True, do_print=False)
    time.sleep(0.2)
    assert in_interval(timer.elapsed_time_last_reset(), (0.2, 0.22))
    time.sleep(0.2)
    time_since_last_call, time_since_last_reset = timer(do_print=False)
    assert in_interval(time_since_last_call, (0.4, 0.42))
    assert in_interval(time_since_last_reset, (0.4, 0.42))
    
def test_cur_time_ns():
    # less than 1 microsecond deviation
    # print(cur_ns_time())
    time_difference = time.time_ns()/1_000_000_000 - cur_ns_time()
    assert abs(time_difference) < 1e-2, f"observed time difference {time_difference}s larger than 1e-2s"
    
    times = []
    for _ in range(100):
        time.sleep(0.01)
        times.append(cur_ns_time())
    assert is_sorted(times)
    assert all(in_interval(x, (0.008, 0.03)) for x in np.diff(times)), f"Max difference {np.diff(times).max()}" # may fail occasionally due to other processes running on the system