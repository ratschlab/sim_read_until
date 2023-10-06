"""
Utilities for the simulator
"""

#%%
from contextlib import contextmanager
import logging
from pathlib import Path
import time
from itertools import count as _count
import numpy as np

from simreaduntil.shared_utils.logging_utils import temp_logging_level


def format_percentage(x):
    """Interprets number in 0, 1 as percentage with two decimal places"""
    eps = 1e-8
    assert -eps <= x <= 1 + eps, f"{x} not in range"
    x = np.clip(x, 0, 1)
    return f"{x * 100:.2f}"

def in_interval(x, interval):
    """
    Check whether a number is in an interval
    
    Args:
        x: number, not a list or vector
        interval: tuple of two numbers
    """
    return interval[0] <= x <= interval[1]

# pylint: disable=invalid-name
_counter = _count()
def new_thread_name(template_str="ont-sim-{}"):
    """
    Helper to generate new thread names
    
    Args:
        template_str: string with one placeholder for the counter
    """
    return template_str.format(next(_counter))

@contextmanager
def set_package_log_level(log_level=logging.INFO):
    """
    Context manager to enable logging for the whole package and ReadFish by setting the levels of the loggers

    This is useful when this module is imported in another script.
    You may also consider calling enable_simulation_logging().__enter__() in a Jupyter notebook.
    """
    # note: child loggers may still have different log levels (e.g. debug which get passed through despite the parent log level)
    with temp_logging_level(logging.getLogger("ru"), log_level):
        with temp_logging_level(logging.getLogger("simreaduntil"), log_level):
            yield
