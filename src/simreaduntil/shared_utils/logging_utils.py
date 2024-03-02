"""
Utilities to use the native logging module

This module is not named logging.py because it would conflict with the logging module
"""

from contextlib import contextmanager
import logging
from pathlib import Path
from typing import Union
import warnings

def setup_logger_simple(name, level=logging.NOTSET):
    """
    Set up logger at module level
    
    Put the following at the top of a file:
    logger = setup_logger_simple(__name__)
    
    NOTSET by default, so it is the level of the first parent which has a level set. If no parent has, it is WARNING.
    
    Args:
        name: use __name__ for the module logger, None for the root logger
        level: logging level (all messages below are filtered)
    
    Returns:
        logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    return logger

module_logger = setup_logger_simple(__name__)
"""module logger"""

def print_logging_levels():
    """Utility function to log current levels of loggers and their handlers"""
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict if (name.startswith("simreaduntil") or name.startswith("ru") or name == "__main__")] + [logging.getLogger()]
    loggers = sorted(loggers, key=lambda x: x.level)
    logger_levels = "\n".join([f"{logger.name}: {logging._levelToName[logger.level]}, propagate={logger.propagate}" for logger in loggers])
    handler_levels = "Handler levels:\n" + "\n".join([f"{logger.name}: " + "; ".join([f"{h.__class__.__name__} with level {logging._levelToName[h.level]}" for h in sorted(logger.handlers, key=lambda x: x.level)]) for logger in loggers if len(logger.handlers) > 0])
    print("Current logging levels:")
    print(logger_levels)
    print(handler_levels)
    module_logger.info("Current logging levels:")
    module_logger.info(logger_levels)
    module_logger.info(handler_levels)

def make_handler_support_end(handler: logging.Handler):
    """
    Make handler support "end" keyword (as in print())
    Use it in conjunction with END_WITH_CARRIAGE_RETURN
    
    Args:
        handler: logging handler (that is responsible for emitting/outputting log records)
    """
    old_emit = handler.emit
    def custom_emit(record):
        # logging.Handler class has a terminator attribute, so we just modify it
        old_terminator = handler.terminator
        handler.terminator = record.__dict__.get("end", "\n")
        old_emit(record)
        handler.terminator = old_terminator
    handler.emit = custom_emit

END_WITH_CARRIAGE_RETURN = {"extra": {"end": "\r"}}
"""use together with make_handler_support_end to print to logger with a carriage return (move to beginning of line, overwriting content)"""

def logging_output_formatter(handler):
    """configures handler to use a specific formatter"""
    formatter = logging.Formatter("%(asctime)s - %(message)s --- %(filename)s:%(lineno)d (%(funcName)s) %(levelname)s ##")
    # "--- vscode://%(pathname)s:%(lineno)d - %(levelname)s"
    make_handler_support_end(handler)
    handler.setFormatter(formatter)
    
_STREAM_HANDLER_ATTR_NAME = "COMPREHENSIVE_STREAM_HANDLER"
def add_comprehensive_stream_handler_to_logger(logger: Union[str, logging.Logger, None]=None, level=logging.NOTSET):
    """
    Add a stream handler to logger that prints detailed information about where log event happened
    
    We only add the handler if it has not already been added to the logger (by checking if it has an attribute).
    It also enables logging of warnings issued to the warnings module to a logger called "py.warnings" at level WARNING.
    
    Args:
        logger: logger, defaults to the root logger
        level: logging level of handler; does not change the logging level of the logger
        
    Returns:
        whether a handler was added
    """
    # add attribute to know that we have added the handler

    if logger is None or isinstance(logger, str):
        logger = logging.getLogger(logger)
    # logger.setLevel(level) # commented because we don't want messages from third-party libraries
    if any(hasattr(h, _STREAM_HANDLER_ATTR_NAME) for h in logger.handlers):
        warnings.warn(f"Logger {logger.name} already has handler attached, not reattaching", stacklevel=2)
        return False

    handler = logging.StreamHandler() # outputs to sys.stderr
    logging_output_formatter(handler)
    handler.setLevel(level)
    setattr(handler, _STREAM_HANDLER_ATTR_NAME, True)

    
    logger.addHandler(handler)
    
    logging.captureWarnings(True)
    return True
    
def remove_comprehensive_stream_handler(logger: Union[str, logging.Logger, None]=None):
    """
    Remove stream handler from logger, counterpart to add_comprehensive_stream_handler_to_logger
    
    Args:
        logger: logger, defaults to the root logger
    """
    if logger is None or isinstance(logger, str):
        logger = logging.getLogger(logger)
    for h in logger.handlers:
        if hasattr(h, _STREAM_HANDLER_ATTR_NAME):
            logger.removeHandler(h)
        
    logging.captureWarnings(False)
    
@contextmanager
def ComprehensiveStreamHandler(logger=None, level=logging.NOTSET):
    """Add a stream handler and restore afterwards again"""
    was_added = add_comprehensive_stream_handler_to_logger(logger, level=level)
    yield
    if was_added:
        remove_comprehensive_stream_handler(logger)
        
def add_filelogging(log_dir, modules):
    """
    Add handlers to log each logger separately to a file
    
    One log file per module is created named as the module name, and one for the root logger.
    The modules that log can be obtained from (typically custom loggers are the last ones in the list):
    logging.getLogger(None).manager.loggerDict.keys()
    
    Args:
        log_dir: directory where to put the log files
        modules: modules to log
    """
    log_dir = Path(log_dir)
    if not log_dir.exists():
        log_dir.mkdir()

    formatter = logging.Formatter("%(asctime)s - %(message)s --- %(filename)s:%(lineno)d (%(funcName)s) %(levelname)s##")
    for module in modules:
        logger = logging.getLogger(module)
        handler = logging.FileHandler(log_dir / (module + ".txt"))
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # root logger
    formatter_root = logging.Formatter("%(asctime)s - %(message)s --- %(filename)s:%(lineno)d (%(funcName)s) %(module)s %(levelname)s##") # adds module
    logger = logging.getLogger(None)
    handler = logging.FileHandler(log_dir / "root.txt")
    handler.setFormatter(formatter_root)
    logger.addHandler(handler)

    def _inform_about_log():
        logger.info(f"""Root log available at: {(log_dir / "root.txt").resolve()}""")
    import atexit
    atexit.register(_inform_about_log)

@contextmanager
def temp_logging_level(logger, level):
    """
    Context manager to set the logging level temporarily and then restore it to the old level
    """
    old_level = logger.level
    logger.setLevel(level)
    yield
    logger.setLevel(old_level)
    