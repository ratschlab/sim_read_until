
import pytest
import logging
import time

from simreaduntil.shared_utils.logging_utils import ComprehensiveStreamHandler, setup_logger_simple, make_handler_support_end, END_WITH_CARRIAGE_RETURN, temp_logging_level

def test_logging_utils():
    logger = setup_logger_simple("test_logging_utils1", level=logging.INFO)
    logger2 = setup_logger_simple("test_logging_utils2", level=logging.INFO)
    
    with ComprehensiveStreamHandler(level=logging.INFO):
        logger.info("hello")
        logger2.info("hello2") # not logged to console because nothing attached to logger2
        
        log_format = 'PREFIX%(message)sSUFFIX'
        formatter = logging.Formatter(log_format)
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        make_handler_support_end(handler)
        logger.addHandler(handler)
        for i in range(5):
            logger.info(f"step {i}", **END_WITH_CARRIAGE_RETURN)
            time.sleep(0.15)

def test_temp_logging_level(caplog):
    logger = logging.getLogger("test_temp_logging_level")
    with temp_logging_level(logger, logging.DEBUG):
        logger.debug("debug message")
        logger.info("info message")
    logger.debug("ignored debug message")
    
    assert caplog.messages == ["debug message", "info message"]
    