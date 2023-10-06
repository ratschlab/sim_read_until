"""
Handle threads, e.g., to see whether they raised an error, or to make a class thread-safe
"""

import contextlib
from functools import wraps
import threading
from threading import RLock
import types

from simreaduntil.shared_utils.logging_utils import setup_logger_simple

logger = setup_logger_simple(__name__)
"""module logger"""

class ThreadWithResultsAndExceptions(threading.Thread):
    """
    Thread that executes a function, gives access to its return value, and saves and logs any exceptions.
    
    This allows to check later whether a thread raised an exception, which otherwise may go unnoticed.
    Calling .join() raises an exception if an exception has occurred in the thread.
    
    """
    def run(self):
        self.exc = None
        try:
            self.ret = self._target(*self._args, **self._kwargs)
        except BaseException as e:
            self.exc = e
            logger.exception("Exception has occurred in thread") # also gives details about exception
            # reraise so that excepthook can take care of it as well, by default just printed (immediately) without aborting the program
            raise e
            
    def join(self, timeout=None):
        """join thread, raising an exception if an exception occurred in the thread"""
        super().join(timeout=timeout)
        self.raise_if_error()
        
    def raise_if_error(self):
        """raise an exception if an exception occurred in the thread (does not wait for thread to finish)"""
        if hasattr(self, "exc") and self.exc is not None:
            logger.error(f"Raising because an exception occurred in thread: {self.exc}")
            raise self.exc

def create_thread_safe_maker():
    """
    Wrapper to apply to a function, the function must first acquire the instance lock self._lock before executing
    """
    def wrapper(f):
        @wraps(f)
        def inner_wrapper(self, *args, **kwargs):
            with self._lock:
                return f(self, *args, **kwargs)
        return inner_wrapper
    return wrapper

class MakeThreadSafe:
    """
    Inheriting from this class makes instance method and instance attribute access thread-safe.
    
    At most one thread at a time can call this class (and it can call recursive methods).
    It can be used as a mixin.
    
    Do not forget to call "super().__init__()".
    
    Modified from https://stackoverflow.com/questions/39145796/locking-a-method-in-python
    and https://www.oreilly.com/library/view/python-cookbook/0596001673/ch06s06.html
    because there, only methods are made thread-safe, not attributes.
    Also added possibility to delete lock for pickling.
    """
    def __init__(self):
        self._lock = RLock()
        
    def __init_subclass__(cls, **kwargs):
        # called as a class method once subclass was initialized, wraps a lock around each instance method
        # use this to avoid having to decorate each instance method individually
        super().__init_subclass__(**kwargs)
        
        make_thread_safe = create_thread_safe_maker()
        for name in cls.__dict__:
            attr = getattr(cls, name)
            if callable(attr) and isinstance(attr, types.FunctionType) and name != "__init__":
                # only add lock to instance methods, not static methods or class methods
                setattr(cls, name, make_thread_safe(attr))
                
    def __getattribute__(self, name):
        # decorate attribute access (including function access, but they are wrapped on top with __init_subclass__)
        
        if name in ["_lock", "__setstate__"]:
            # support pickling
            return object.__getattribute__(self, name)
        
        # with self._lock: # triggers __getattribute__ again
        with object.__getattribute__(self, "_lock"):
            # print(f"Thread {threading.get_ident()}: Acquired lock to get attribute {name}")
            return object.__getattribute__(self, name)
                
    def _delete_lock(self):
        """Remove lock, useful for pickling"""
        self._lock = contextlib.nullcontext() # context manager that does nothing
        
    # support pickling
    def __getstate__(self):
        state = self.__dict__.copy()
        del state["_lock"]
        return state
    def __setstate__(self, state):
        self._lock = RLock()
        self.__dict__.update(state)
        