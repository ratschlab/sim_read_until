"""
Functions for profiling and logging memory usage
"""

import tracemalloc

from simreaduntil.shared_utils.utils import setup_logger_simple

from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger

logger = setup_logger_simple(__name__)
"""module logger"""

########## Memory usage (full Python program since tracing was started): peak + current size ##########

def set_profiling(enable):
    """Enables or disables profiling"""
    if enable:
        tracemalloc.start()
    else:
        tracemalloc.stop()
        
class ProfilingEnabled():
    """
    Make sure that profiling is enabled. If it was not already, it will be disabled after the "with" block.
    
    Does not do anything if profiling is already running (e.g. recursive calls) since not possible to restore properly
    """
    def __enter__(self):
        self.was_running = tracemalloc.is_running()
        if not self.was_running:
            tracemalloc.start()
    def __exit__(self):
        if not self.was_running:
            tracemalloc.stop()
        
# see https://stackoverflow.com/questions/1094841/get-human-readable-version-of-file-size
def _sizeof_fmt(num, suffix="B"):
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:8.4f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.5f}Yi{suffix}"

# program should have been started with PYTHONTRACEMALLOC=1 or use tracemalloc.start()/.stop()
# it logs it under the function that calls this function (when reporting function name etc.)
def log_memory_usage(prefix=""):
    """
    Log memory usage (since tracemalloc was started)
    
    Args:
        prefix: to prefix the log message with, for easier identification in log file
    """
    if tracemalloc.is_tracing():
        # unit is bytes: see https://github.com/python/cpython/blob/95e271b2266b8f2e7b60ede86ccf3ede4a7f83eb/Lib/test/test_tracemalloc.py#L233
        current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        logger.debug(f"{prefix}Memory: Current size: {_sizeof_fmt(current_bytes)}, peak size: {_sizeof_fmt(peak_bytes)}", stacklevel=2)
        return current_bytes, peak_bytes
    return None


########## approximate memory usage by summing up object sizes, avoids double counting the same object (using id()) ##########
# may not take into account everything of an object (just an int pointer in the worst case, rather than the actual object size)

from sys import getsizeof, stderr
from itertools import chain
from collections import deque, defaultdict

# print approximate size of an object, including its contents (not just its pointers), in bytes
# see https://code.activestate.com/recipes/577504/
def total_approx_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}
                    
    Examples:
        d = dict(a=1, b=2, c=3, d=[4,5,6,7], e='a string of chars')
        total_size(d, verbose=True)

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    defaultdict: dict_handler, # added
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            logger.debug(s, type(o), repr(o), file=stderr)

        handled = False
        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                handled = True
                break
        if not handled:
            # this should deal with most user-defined classes
            # assumes no __slots__ used
            if hasattr(o, '__dict__'):
                s += sizeof(o.__dict__)
        return s

    return sizeof(o)

def log_object_sizes(*args, **kwargs):
    """
    Estimate the size of an object and log it to the logger
    
    Args:
        args: unnamed objects, will be labeled by index in output
        kwargs: named objects, will be labeled by key in output
        
    Example:
        log_object_sizes(a=1, b=3)
    """
    all_args = {f"arg{i}": arg for (i, arg) in enumerate(args)}
    all_args.update(kwargs)
    logger.debug(f"Object sizes: " + ", ".join("{}: {}".format(name, _sizeof_fmt(total_approx_size(o))) for (name, o) in all_args.items()), stacklevel=2)


if __name__ == "__main__":
    add_comprehensive_stream_handler_to_logger()
    
    aa = {"aa": 3, "bb": 5}
    
    log_object_sizes(a=1, b=3)
    log_object_sizes(33, aa=aa, a=1)