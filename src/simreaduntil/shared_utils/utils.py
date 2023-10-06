"""
General utility functions
"""

import copy
import functools
import gzip
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, Iterable, List
import dill
import tqdm

from simreaduntil.shared_utils.logging_utils import setup_logger_simple

        
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