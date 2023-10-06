"""
Tee stdout and stderr to a file
"""

from io import TextIOBase
import logging
import sys
import threading

class TeeStdouterr:
    """
    Tee stdout and stderr
    
    We use the logging module to ensure thread-safety when writing to it. Do not use this logger or its children for anything else.
    
    Unlike common implementations seen on the Internet, this implementation correctly logs stdout and stderr to the respective 
    streams (instead of merging them into one stream).
    
    Args:
        filename (str): file to tee to
        mode (str): file mode
        propagate (bool): whether to propagate to root logger; if set to True, every write to stdout and stdin (e.g. prints) can be intercepted with the usual logging additionally
    """
    def __init__(self, filename, mode="w", propagate=False):
        logger = logging.getLogger(TeeStdouterr.LOGGER_NAME)
        logger.setLevel(logging.INFO)
        logger.propagate = propagate
        file_handler = logging.FileHandler(filename, mode=mode)
        file_handler.terminator = ""
        logger.addHandler(file_handler)
        
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.terminator = ""
        logger.getChild("stdout").addHandler(stream_handler)
        stream_handler = logging.StreamHandler(sys.stderr)
        stream_handler.terminator = ""
        logger.getChild("stderr").addHandler(stream_handler)
        
        self._stdout = None
        self._stderr = None
    
    LOGGER_NAME = "stdouterr"
    """Name of logger (do not interact with this logger directly)"""
    
    def redirect(self):
        """Start redirection"""
        assert self._stdout is None and self._stderr is None, "already started"
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self.Writer("stdout")
        sys.stderr = self.Writer("stderr")
    
    def restore(self):
        """
        Restore stdout and stderr
        """
        assert self._stdout is not None and self._stderr is not None, "not started"
        sys.stdout = self._stdout
        sys.stderr = self._stderr
        self._stdout = None
        self._stderr = None
    
    def is_redirecting(self):
        """
        Whether teeing is currently active
        """
        assert (self._stdout is not None) == (self._stderr is not None)
        return self._stdout is not None
    
    class Writer(TextIOBase):
        """
        Writer class that writes to the corresponding child logger which gets propagated to its parent
        
        Args:
            what: child logger to write to
        """
        def __init__(self, what):
            self._file_logger = logging.getLogger(TeeStdouterr.LOGGER_NAME)
            self._stream_logger = self._file_logger.getChild(what)
    
        # write interface
        
        def write(self, data):
            # logging.Logger is thread-safe
            self._stream_logger.info(data)
        
        def flush(self):
            self._file_logger.handlers[0].flush()
            self._stream_logger.handlers[0].flush()
        
        def close(self):
            logger = self._file_logger
            if logger is not None:
                logger.handlers[0].close()
                logger.removeHandler(logger.handlers[0])
                self._file_logger = None


if __name__ == "__main__":
    # todo2: move
    
    import time

    # should not appear in the file
    print("hello")
    print("hello2", file=sys.stderr)

    out_filename = "test.txt"
    tee_object = TeeStdouterr(out_filename)
    tee_object.redirect()

    def print_to_stdout():
        for i in range(1000):
            # print(i)
            sys.stdout.write(f"{i}\n")
            time.sleep(1e-5)
            
    def print_to_stderr():
        for i in range(2000, 2000+1000):
            # todo2
            # print(i, file=sys.stderr)
            sys.stderr.write(f"{i}\n")
            time.sleep(1e-5)

    t1 = threading.Thread(target=print_to_stdout)
    t2 = threading.Thread(target=print_to_stderr)
    t1.start()
    t2.start()
    t1.join()
    t2.join()

    sys.stdout.flush()
    sys.stderr.flush()
    # tee_object.restore()

    # read file and check all numbers appear
    with open(out_filename, "r") as f:
        lines = [l.strip() for l in f.readlines()]
        expected_lines = set(map(str, range(1000))).union(set(map(str, range(2000, 2000+1000))))
        assert len(lines) == len(expected_lines)
        assert set(lines) == expected_lines