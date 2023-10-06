import pytest

import sys
import time

from simreaduntil.shared_utils.thread_helpers import ThreadWithResultsAndExceptions, MakeThreadSafe
from simreaduntil.shared_utils.utils import dill_dump, dill_load

def test_thread_with_noexceptions():
    def f1():
        print(1)
        time.sleep(0.2)
        print(2)
        return "hello"
    
    x = ThreadWithResultsAndExceptions(target=f1, daemon=True)
    x.start()
    x.join()
    assert x.ret == "hello"
    assert x.exc == None
    
    x = ThreadWithResultsAndExceptions(target=f1, daemon=True)
    x.start()
    time.sleep(0.5)
    x.join(timeout=0)
    assert x.ret == "hello"
    
# silence pytest.PytestUnhandledThreadExceptionWarning since thread raises an exception
@pytest.mark.filterwarnings("ignore::pytest.PytestUnhandledThreadExceptionWarning")
def test_thread_with_exceptions():
    # this example should print 4 exception traces (for logger.exception, raise e, logger.error, raise self.exc)
    def f_with_error():
        print(1)
        time.sleep(0.2)
        raise ValueError("some error") # note: vscode goes into debugger, so run from pytest sidebar panel
        return "hello" # never reached
    
    x = ThreadWithResultsAndExceptions(target=f_with_error, daemon=True)
    x.start()
    sys.stdout.flush()
    x.join(timeout=0.05) # not raising exception yet
    assert x.exc is None
    
    with pytest.raises(ValueError, match="some error"):
        x.join()
    assert x.exc is not None
    assert not hasattr(x, "ret")

def test_MakeThreadSafeSimple(tmp_path):
    class SimpleExample(MakeThreadSafe):
        def __init__(self):
            super().__init__()
            self.a1 = 2
            
        def hello(self):
            print("hello")
            
    obj = SimpleExample()

    assert obj.a1 == 2
    obj.hello()
    obj._delete_lock()
    
    # test pickling
    obj = SimpleExample()
    obj.a1 = 3
    pickle_filename = tmp_path / "test.dill"
    dill_dump(obj, pickle_filename)
    obj = dill_load(pickle_filename)
    obj.hello()
    assert obj.a1 == 3
    
def test_MakeThreadSafeConcurrency():
    # test that locking properly works
    
    class ExampleThreadSafeClass(MakeThreadSafe):
        def __init__(self) -> None:
            super().__init__()
            self.i = 0
            
        def hello(self):
            print("Starting now"); sys.stdout.flush()
            v = self.i
            time.sleep(0.5)
            self.i = v + 1
            print(f"Hello {self.i}"); sys.stdout.flush()
            return self.i
            
    syncObj = ExampleThreadSafeClass()
    t1 = ThreadWithResultsAndExceptions(target=syncObj.hello)
    t2 = ThreadWithResultsAndExceptions(target=syncObj.hello)
    
    t1.start()
    time.sleep(0.02) # so that t1 acquires the lock first
    t2.start()
    
    t1.join()
    t2.join()
    
    assert t1.ret == 1
    assert t2.ret == 2
    assert syncObj.i == 2
    
    # remove lock
    syncObj._delete_lock()
    t1 = ThreadWithResultsAndExceptions(target=syncObj.hello)
    t2 = ThreadWithResultsAndExceptions(target=syncObj.hello)
    t1.start()
    t2.start()
    t1.join()
    t2.join()
    assert t1.ret == 3
    assert t2.ret == 3
    assert syncObj.i == 3