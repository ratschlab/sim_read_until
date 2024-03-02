"""
General utility functions
"""

from contextlib import contextmanager
import contextlib
import copy
import functools
import gzip
import os
from pathlib import Path
import queue
import shutil
import subprocess
import sys
import threading
import time
from typing import Any, Dict, Iterable, List
import dill
import tqdm

from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.timing import cur_ns_time

        
logger = setup_logger_simple(__name__)
"""module logger"""

###### general utils

def is_sorted(arr: List[Any]):
    """Whether a list is sorted"""
    return all(x <= y for (x, y) in zip(arr[:-1], arr[1:]))

def subset_dict(d: Dict, keys: Iterable[Any]):
    """Subset dict to some keys, all keys must be present"""
    return {k: d[k] for k in keys}

def get_some_value_from_dict(d: Dict):
    """
    Get some value from a dict
    
    Raises:
        StopIteration if dict is empty
    """
    return next(iter(d.values()))

def convert_to_vscode_remote_path(path):
    """Convert path on remote to a path that VSCode can open"""
    return "vscode-remote://" + str(path)

###### directory methods

def delete_dir_if_exists(dir, ask=False):
    """
    Delete directory if it exists, return whether it existed
    
    Args:
        dir: directory to delete
        ask: whether to ask the user before deleting if the directory exists; if user rejects, the directory is not deleted
        
    Returns:
        whether the directory was deleted
    """
    dir = Path(dir)
    if dir.exists():
        if (not ask) or input(f"Delete directory '{dir.resolve()}'? y/n: ").strip().lower() == "y":
            shutil.rmtree(dir)
            return True
    return False

def is_empty_dir(path):
    """Whether the directory is empty"""
    for _ in Path(path).iterdir():
        return False
    return True

def num_lines_in_file(filename):
    """
    Counts empty lines, stripping last line if empty
    """
    open_method = functools.partial(gzip.open, mode="rb") if str(filename).endswith(".gz") else functools.partial(open, mode="r")
    with open_method(filename) as f:
        return sum([1 for _ in f])
    
def force_eval_generator(gen):
    """
    Evaluate a generator to its first yield.
    
    Generator functions (using yield in the function) are only run when the first element is requested.
    This may come with some overhead for the first element, e.g. to open a FASTA file. This function
    avoids it by evaluating the generator to its first yield.
    
    Args:
        gen: a generator
        
    Returns:
        new generator object (do not use the old one which has one element missing)
    """
    import inspect
    from more_itertools import peekable
    assert inspect.isgenerator(gen), "must be a generator"
    gen = peekable(gen)
    gen.peek()
    return gen

def force_eval_generator_function(f):
    """
    Force evaluate a generator function until its first yield, decorate a function with @force_eval_generator_function
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return force_eval_generator(f(*args, **kwargs))
    return wrapper

def dill_dump(obj, filename):
    """
    Dumps an object to a file using dill.
    
    Dill has the advantage over pickle that it can serialize lambdas and other functions.
    """
    with open(filename, "wb") as f:
        dill.dump(obj, f)
        
def dill_load(filename):
    """
    Load an object from a file using dill.
    """
    with open(filename, "rb") as f:
        return dill.load(f)
    
def get_file_content(filename):
    """return the file content as a string"""
    with open(filename) as f:
        return f.read()
    
def is_empty_file(filename):
    """whether the file is empty"""
    return os.stat(filename).st_size == 0

# verbose=False -> returns stdout (and stderr merged into it)
def print_cmd_and_run(cmd, verbose=True, dry=False, **kwargs):
    """
    Print the command and run it, raising an error if it fails
    
    Args:
        cmd: command to run
        verbose: whether to print the command
        dry: whether to actually run the command, will print the command
        kwargs: to pass to subprocess.run
        
    Returns:
        stdout if verbose=False, otherwise None
    """
    if dry:
        logger.warning(f"Dry run, so not executing the command:\n{cmd}", stacklevel=2)
    else:
        if verbose:
            logger.info(f"Executing the following command:\n{cmd}", stacklevel=2)
        try:
            p = subprocess.run(cmd, shell=True, check=True, capture_output=not verbose, **kwargs)
        except subprocess.CalledProcessError as e:
            if not verbose:
                # stdout, stderr only available when capture_output is False
                logger.error(f"Process failed with stdout:\n{e.stdout.decode()}\nand stderr:\n{e.stderr.decode()}")
            raise e
        
        if not verbose:
            return p.stdout.decode()
        # else: cannot both print to the console and capture output easily

def tqdm_if_verbose(iter, *args, verbose=False, **kwargs):
    """
    Return tqdm if verbose, otherwise just the iterator
    """
    if verbose:
        return tqdm.tqdm(iter, *args, **kwargs)
    else:
        return iter

def tqdm_with_name(it, *args, **kwargs):
    """
    Use this instead of tqdm.tqdm
    
    Arguments:
        it: iterator returning (name, elem) rather than just elem, name will be displayed in progress bar 
            when it is available and when the next element is requested
            
        args, kwargs: to pass to tqdm.tqdm
    """
    # it has returns name and elem (once elem is returned)
    progress_bar = tqdm.tqdm(it, *args, **kwargs)
    for (name, elem) in progress_bar:
        progress_bar.set_description(f"computed {name}, processing")
        yield elem
        progress_bar.set_description(f"finished {name}, getting next")
        
def print_args(args, logger=None, exclude=None):
    """
    Print and format all arguments from the command line
    
    Copied from ReadFish
    """
    if exclude is None:
        exclude = []
    dirs = dir(args)
    m = max([len(a) for a in dirs if a[0] != "_"])
    for attr in dirs:
        if attr[0] != "_" and attr not in exclude and attr.lower() == attr:
            record = "{a}={b}".format(a=attr, m=m, b=getattr(args, attr))
            if logger is not None:
                logger.info(record)
            else:
                print(record)

def _cycle_list_deep(lst):
    """
    Cycle through a list, returning the next element each time

    Args:
        lst: list to cycle through

    Yields:
        next element
    """
    # itertools cycle may return the same element several times, here we copy it instead
    import copy
    while True:
        for x in lst:
            yield copy.deepcopy(x)


@contextmanager
def set_signal_handler(signal_type, handler):
    """
    Set a signal handler temporarily

    This function can only be called from the main thread of the main interpreter.

    This is better than using KeyboardInterrupt because this immediately breaks out of the code, possibly leaving the code in an inconsistent 
    state, e.g. for the simulator updating channels when Ctrl-C is pressed causes it to immediately stop. It is better to set a flag to stop
    it. Since the signal handler acquires the Pytho GIL, make sure the code in the signal handler can run to completion because all other threads are 
    stopped as well due to the Python GIL, so no locked mutexes should be required.

    See the tests for an example.

    Args:
        signal_type: e.g. signal.SIGINT for KeyboardInterrupt
        handler: function to handle the signal taking two arguments, should not mess with any state
    """
    import signal
    old_handler = signal.getsignal(signal_type)
    signal.signal(signal_type, handler)
    yield
    signal.signal(signal_type, old_handler)


@contextmanager
def tee_stdouterr_to_file(filename_base, mode="a"):
    """
    Try to use a file handler with the logger rather than this function!

    Tee the entire output of a Python program to a file and the console.

    # assign to an object since otherwise it will get destroyed and the file handler will become invalid!
    obj = tee_stdouterr_to_file("test11", mode="w")
    obj.__enter__()
    # or use with a "with" statement

    Args:
        filename_base: base of the filename to write to, ".out" and ".err" will be appended, directory containing the file must exist
        mode: mode to open file in, "a" for append, "w" for overwrite
    """

    """
    File object that writes to multiple file objects at once, e.g. for teeing
    """
    class MultipleFilesWriter:
        def __init__(self, *files):
            self.files = files

        def write(self, text):
            for file in self.files:
                file.write(text)

        def flush(self):
            for file in self.files:
                file.flush()

    # to debug this function in case of exception, remove the redirection of stdout
    with open(str(filename_base) + ".out", mode=mode) as out_file:
        with open(str(filename_base) + ".err", mode=mode) as err_file:
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            with contextlib.redirect_stdout(MultipleFilesWriter(old_stdout, out_file)):
                with contextlib.redirect_stderr(MultipleFilesWriter(old_stderr, err_file)):
                    yield


"""
Class storing a value

Useful for passing values by reference to a function, since int, float etc are immutable, so modifications to them are not visible outside the function
"""
class MutableValue:
    def __init__(self, value=None):
        self.value = value
# def __repr__(self):
#     return f"MutableValue({self.value})"
# def __str__(self):
#     return f"MutableValue({self.value})"
# def __eq__(self, other):
#     return self.value == other.value
# def __hash__(self):
#     return hash(self.value)

def record_gen_waiting_time(gen, waiting_time: MutableValue):
    """
    Record how much time is spent waiting for new elements in the generator, only counting the times when requesting an element and until getting it

    The function also works if the generator is stopped early.
    
    Args:
        gen: generator to wrap
        waiting_time: MutableValue to store the total waiting time in seconds, only includes the time until the generator is destroyed
    
    Yields: values from gen
    """
    elapsed_time = 0
    try:
        t_before = cur_ns_time()
        for x in gen:
            elapsed_time += cur_ns_time() - t_before
            yield x
            t_before = cur_ns_time()
    finally:
        # to make sure this is executed even if the generator is stopped early
        waiting_time.value = elapsed_time
        
def record_gen_fcn_waiting_time(gen_fcn, gen, waiting_time: MutableValue):
    """
    Record how much time the function gen_fcn(gen) takes processing elements coming from generator gen, excluding the time waiting for elements from gen
    
    Example:
    See tests
    
    Args:
        gen_fcn: function taking a generator and returning a generator
        gen: generator to wrap
        waiting_time: MutableValue to store the total waiting time in seconds of gen_fcn itself only, only includes the time until the generator is destroyed
    """
    waiting_time.value = 0
    time_inner_gen = MutableValue()
    try:
        yield from record_gen_waiting_time(gen_fcn(record_gen_waiting_time(gen, time_inner_gen)), waiting_time)
    finally:
        waiting_time.value -= time_inner_gen.value
        

"""
Queue with interruptible get and put methods when it is stopped, returned a QueueStoppedException
"""
class StoppableQueue(queue.Queue):
    class QueueStoppedException(Exception):
        pass
    # workaround since otherwise Empty not found
    Empty = queue.Empty
    Full = queue.Full
    
    def __init__(self, maxsize=0):
        super().__init__(maxsize)
        # Create an event that can be used to signal the queue to stop its operations.
        self.stop_event = threading.Event()
        self.POLL_INTERVAL = 0.1

    """
    Put an item on the queue
    
    raises an QueueStoppedException if the stop event is set
    
    Args:
        item: the item to put on the queue
        block: if True, repeatedly try to put item (repeating if queue is full); if False, immediately raise queue.Full if full
        timeout: only applies when block=True
    """
    def put(self, item, block=True, timeout=None):
        start_time = time.time()
        while not self.stop_event.is_set():
            try:
                super().put(item, block, timeout=self.POLL_INTERVAL)
                return
            except queue.Full:
                if (not block) or (timeout is not None and time.time() - start_time > timeout):
                    raise # queue full exception
        if self.stop_event.is_set():
            raise self.QueueStoppedException()

    """
    Get an item from the queue
    
    raises an QueueStoppedException if the stop event is set
    
    Args:
        block: if True, repeatedly try to get item (repeating if queue is empty); if False, immediately raise queue.Empty if empty
        timeout: only applies when block=True
        
    Returns:
        item from the queue, or exception if None (queue empty or stop event set or timeout expired)
    """
    def get(self, block=True, timeout=None):
        start_time = time.time()
        while not self.stop_event.is_set():
            try:
                item = super().get(block, timeout=self.POLL_INTERVAL)
                return item
            except queue.Empty:
                if (not block) or (timeout is not None and time.time() - start_time > timeout):
                    raise  # queue empty exception
        if self.stop_event.is_set():
            raise self.QueueStoppedException()

    """
    Convenience method to set the stop event.
    """
    def stop(self):
        self.stop_event.set()