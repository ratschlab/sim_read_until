"""
Helper functions for debugging
"""

import os

from simreaduntil.shared_utils.logging_utils import setup_logger_simple

logger = setup_logger_simple(__name__)
"""module logger"""

def is_test_mode():
    """
    Detect if test mode is active
    Scripts can be tested before the main run by setting the environment variable 
    export TEST_MODE=1
    If test mode is active, it executes something else and automatically warns about debugging.
    
    Returns:
        whether the program is run in test mode.
        
    Example:
        if is_test_mode():
            print("In test mode")
        else:
            print("In normal mode")
            
    """
    val = os.environ.get("TEST_MODE", "0").strip() == "1"
    if val:
        warn_debugging()
    return val

__WARN_DEBUGGING_REGISTERED = False # to only warn once
def warn_debugging():
    """
    Warn about debugging if not already done so (using print() and logger)
    Use it in conjunction with .is_test_mode().
    """
    def helper():
        # also print, in case logging is disabled
        [print("#"*80) for _ in range(5)]
        print("Running in debug mode")
        [print("#"*80) for _ in range(5)]

        [logger.info("#"*80) for _ in range(5)]
        logger.info("Running in debug mode", stacklevel=2)
        [logger.info("#"*80) for _ in range(5)]
    helper()
    # also print when program terminates, but only once
    global __WARN_DEBUGGING_REGISTERED
    if not __WARN_DEBUGGING_REGISTERED:
        import atexit
        atexit.register(helper)
        __WARN_DEBUGGING_REGISTERED = True
        
if __name__ == "__main__":
    
    prev_val = os.environ.get("TEST_MODE", None)
    
    os.environ["TEST_MODE"] = "0"
    assert not is_test_mode()
    del os.environ["TEST_MODE"]
    assert not is_test_mode()
    
    os.environ["TEST_MODE"] = "1"
    assert is_test_mode()
    
    if prev_val is None:
        del os.environ["TEST_MODE"]
    else:
        os.environ["TEST_MODE"] = prev_val