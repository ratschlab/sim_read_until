"""
Functions related to timing
"""

import time

class ElapsedTimer:
    """
    include_sleep=True -> time.sleep() not included for example
    """
    def __init__(self, include_sleep=True):
        if include_sleep:
            self.timer = time.perf_counter_ns
        else:
            self.timer = time.process_time_ns
        self.last_time = self.timer() # time in nanoseconds, only relative time difference is correct

    # reset: whether to reset elapsed time to zero
    def __call__(self, reset=True, print_time=False):
        cur_time = self.timer()
        elapsed_time = cur_time - self.last_time
        if reset:
            self.last_time = cur_time

        elapsed_time = elapsed_time / 1_000_000_000
        if print_time:
            print(f"{elapsed_time} s elapsed")
        return elapsed_time

class StepAndTotalElapsedTimer:
    """
    Reports incremental (since last call) and total elapsed time (since last reset)
    
    Args:
        include_sleep: whether to include time.sleep() in the elapsed time
    """
    def __init__(self, include_sleep=True):
        self.timer = time.perf_counter_ns if include_sleep else time.process_time_ns
        t = self.timer()
        self.time_last_reset = t
        self.time_last_call = t

    def __call__(self, msg_prefix="", reset=False, do_print=True):
        cur_time = self.timer()
        elapsed_time_last_reset = (cur_time - self.time_last_reset) / 1_000_000_000
        elapsed_time_last_call = (cur_time - self.time_last_call) / 1_000_000_000

        if do_print:
            print(f"{msg_prefix}{elapsed_time_last_call} s (total {elapsed_time_last_reset} s)")
        cur_time = self.timer()
        if reset:
            self.time_last_reset = cur_time
        self.time_last_call = cur_time

        return (elapsed_time_last_call, elapsed_time_last_reset)

    def elapsed_time_last_reset(self):
        """Returns time since last reset without registering it as a call (unlike .())"""
        return (self.timer() - self.time_last_reset) / 1_000_000_000


_time_offset = time.time_ns() - time.perf_counter_ns() # time when this module is loaded, minus offset
def cur_ns_time():
    """
    Get monotonic current time with nanosecond precision
    
    Notes:
    time.time_ns() is not monotonic, so when waking up the machine again, the time might decrease by 1 second or so, so it is not monotonic
    monotonic_ns() and perf_counter_ns() are monotonic, but have an undefined reference point, so only the difference is valid; perf_counter_ns has higher precision    
    
    Returns:
        current time in seconds, slightly shifted if time.time_ns() is not monotonic
    """
    # t0 + (t1_perf - t0_perf), _time_offset = t0 - t0_perf
    return (_time_offset + time.perf_counter_ns())/1_000_000_000