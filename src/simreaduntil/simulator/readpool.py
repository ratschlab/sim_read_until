"""
Read pool that returns reads when requested, e.g., from a generator or a file
"""

import contextlib
from pathlib import Path
from queue import Empty
import queue
import threading
from typing import Optional, Dict, Any, Generator, Tuple
import numpy as np
import pysam
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.thread_helpers import ThreadWithResultsAndExceptions

from simreaduntil.shared_utils.utils import StoppableQueue, force_eval_generator_function, get_some_value_from_dict, is_empty_file

logger = setup_logger_simple(__name__)
"""module logger"""

@force_eval_generator_function
def reads_from_file_gen(fasta_file, shuffle_rand_state: Optional[np.random.Generator]=None):
    """
    Generator returning reads from file

    Also supports .gz files
    
    Yields:
        (id, seq)
    """
    if is_empty_file(fasta_file):
        return
    
    with pysam.FastaFile(fasta_file) as fasta:
        ref_names = fasta.references
        if shuffle_rand_state is not None:
            shuffle_rand_state.shuffle(ref_names) # in-place
        for id in ref_names:
            yield (id, fasta.fetch(id))
    
# todo9: put into readpool
class NoReadLeftException(Exception):
    """
    When no read is left in the read pool
    """
    pass

class ReadPool:
    """
    Read pool from which reads can be obtained
    
    Do not forget to call finish() when done.
    
    Args:
        reads_per_channel: whether reads are channel-specific or any read can be assigned to any channel
    """
    def __init__(self, reads_per_channel):
        self.lock = threading.Lock()
        self.definitely_empty = False # whether the readpool is definitely empty (if False, we don't know)
        self.nb_reads_returned = 0
        self.reads_per_channel = reads_per_channel
        
    def get_new_read(self, channel=None) -> Tuple[str, Any]:
        """
        Get new read (thread-safe)
        
        Once the function returns NoReadLeftException for a channel, it will always return NoReadLeftException for this channel (also for channel=None).
        
        Args:
            channel: channel for which to get read
        """
        with self.lock:
            res = self._get_new_read(channel=channel) if self.reads_per_channel else self._get_new_read()
            self.nb_reads_returned += 1
            return res
    
    def _get_new_read(self, channel=None) -> Tuple[str, Any]:
        """
        Get new read (not thread-safe), overwrite in subclasses
        
        Args:
            channel: channel for which to get read, only provided if self.reads_per_channel is True
            
        Returns:
            tuple (read, read signal)
        """
        raise NotImplementedError()
    
    def finish(self):
        """
        Stop the read pool
        
        For example, if it is threaded, stop the thread
        """
        pass

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.finish()
        
class ReadPoolFromIterable(ReadPool):
    """
    Read pool that requests reads from generator
    """
    def __init__(self, reads_iterable):
        super().__init__(reads_per_channel=False)
        
        self.reads_iterable = reads_iterable
        
    def _get_new_read(self) -> Tuple[str, Any]:
        # note: generators are not thread-safe !!
        try:
            return next(self.reads_iterable)
        except StopIteration as e:
            self.definitely_empty = True
            raise NoReadLeftException from e # exception chaining
        
    # support pickling
    def __getstate__(self):
        # delete generator since a generator cannot be pickled/copied easily
        obj = self.__dict__.copy() # shallow copy, so not copying possible generator
        if isinstance(self.reads_iterable, Generator):
            # cannot pickle a generator
            obj["reads_iterable"] = None
        return obj
        
class ReadPoolFromIterablePerChannel(ReadPool):
    """
    Read pool that requests reads from channel-specific generator
    """
    def __init__(self, reads_iterable_per_channel: Dict[Any, Generator[str, None, None]]):
        super().__init__(reads_per_channel=True)
        
        self.reads_iterable_per_channel = reads_iterable_per_channel
        
    def _get_new_read(self, channel) -> Tuple[str, Any]:
        try:
            return next(self.reads_iterable_per_channel[channel])
        except StopIteration as e:
            raise NoReadLeftException from e # exception chaining
        
    def __repr__(self):
        return f"ReadPoolFromIterable({list(self.reads_iterable_per_channel.keys())})"
        
    # support pickling
    def __getstate__(self):
        # delete generator since a generator cannot be pickled/copied easily
        obj = self.__dict__.copy() # shallow copy, so not copying possible generator
        if len(self.reads_iterable_per_channel) > 0 and isinstance(get_some_value_from_dict(self.reads_iterable_per_channel), Generator):
            # cannot pickle a generator
            obj["reads_iterable_per_channel"] = None
        return obj
        
class ReadPoolFromFile(ReadPoolFromIterable):
    """
    Read pool that reads from a file or directory
    """
    def __init__(self, reads_file_or_dir, shuffle_rand_state: Optional[np.random.Generator]=None):
        reads_file_or_dir = Path(reads_file_or_dir)
        def read_gen():
            for filename in (reads_file_or_dir.glob("**/*.fasta") if reads_file_or_dir.is_dir() else [reads_file_or_dir]):
                    logger.info(f"Starting to read file '{filename}'")
                    yield from reads_from_file_gen(filename, shuffle_rand_state=shuffle_rand_state)
        super().__init__(read_gen())
        
        self.shuffled = shuffle_rand_state is not None
        self.reads_file_or_dir = reads_file_or_dir
        
    @staticmethod
    def can_handle(file: Path) -> bool:
        """
        Check if the read pool can open the file/directory
        
        Args:
            file: file or directory
        """
        file = Path(file)
        return (
            (file.is_dir() and any(file.glob("**/*.fasta"))) or
            (file.is_file() and (file.suffix == ".fasta" or file.suffix == ".fasta.gz"))
        )
        
    def __repr__(self):
        return f"ReadPool(file = {self.reads_file_or_dir}, shuffled = {self.shuffled})"
    
class ThreadedReadPoolWrapper(ReadPool):
    """
    Threaded ReadPool that wraps another ReadPool and reads from it in another thread using a queue

    Note: Using a rng with ThreadedPoolWrapper is not thread-safe if rng is accessed from multiple threads
    """
    
    def __init__(self, read_pool: ReadPool, queue_size: int):
        super().__init__(reads_per_channel=read_pool.reads_per_channel)
        self._read_pool = read_pool
        assert queue_size > 0 # otherwise will read all reads at once
        self._reads_queue = StoppableQueue(queue_size)
        self._reader_thread = ThreadWithResultsAndExceptions(target=self._fill_queue, name="ThreadedReadPoolWrapper")
        self._reader_thread.start()
        self.definitely_empty = False
        
    def can_handle(self, *args, **kwargs) -> bool:
        """
        Check if the read pool can open the file/directory
        
        Params:
            args, kwargs: passed to wrapped read pool
        """
        return self._read_pool.can_handle(*args, **kwargs)
    
    def _fill_queue(self):
        """Keep the queue filled"""
        try:
            while True:
                read = self._read_pool.get_new_read()
                self._reads_queue.put(read)
        except (StoppableQueue.QueueStoppedException, NoReadLeftException):
            pass
        
        try:
            self._reads_queue.put(None) # allows get_new_read() to detect it will never return a read again
        except StoppableQueue.QueueStoppedException:
            # was terminated in between
            pass
        
        logger.info("Finished read queue filler thread")
        
    def __repr__(self) -> str:
        return f"ThreadedReadPool({self._read_pool}, queue_size: {self._reads_queue.maxsize})"
        
    # protected by lock/mutex
    def _get_new_read(self) -> Tuple[str, Any]:
        if self.definitely_empty:
            raise NoReadLeftException
        
        try:
            read = self._reads_queue.get()
        except StoppableQueue.QueueStoppedException:
            read = None
            
        if read is None:
            self.definitely_empty = True
            raise NoReadLeftException
        return read
    
    def finish(self):
        self._reads_queue.stop()
        self._reader_thread.join()
        self._reader_thread.raise_if_error()
        self._read_pool.finish()
        
def get_slow5_reads_gen(filename):
    """generator returning reads in a slow5 file"""
    import pyslow5
    with contextlib.closing(pyslow5.Open(str(filename), "r")) as fh:
        for read in fh.seq_reads():
            yield (read["read_id"], read["signal"])

# """
# Read Slow5 files in another thread and put the read data into a queue.

# # todo: remove read_buffer, rather use ThreadedReadPoolWrapper
# # todo: subclass ReadPool
# # todo: overwrite _get_new_read
# # todo: number of threads for reading slow5

# Args:
#     s5_dir: directory containing slow5 files
#     read_buffer: number of reads to buffer in queue
# """
# class Slow5ReadPool(ReadPool):
#     def __init__(self, s5_dir, read_buffer) -> None:
#         super().__init__(reads_per_channel=False)
        
#         raise ("NotYetImplementedError")
        
#         import pyslow5 # todo: add as dependency
#         from pathlib import Path
#         from queue import Queue
#         import threading
        
#         s5_dir = Path(s5_dir)
#         assert s5_dir.is_dir()
#         self.s5_files = list(s5_dir.glob("*.[sb]low5"))
#         self.queue = Queue(read_buffer) # threadsafe
#         self.cur_file_idx = -1
#         self._queue_filler_thread = threading.Thread(target=self._fill_queue)
#         self._queue_filler_thread.start()
#         # self._queue_filler_thread = multiprocessing.Process(target=self._fill_queue)
        
#     """
#     Check if the read pool can open the file/directory
#     """
#     @staticmethod
#     def can_handle(dir: Path) -> bool:
#         dir = Path(dir)
#         return dir.is_dir() and any(dir.glob("**/*.[sb]low5"))
    
#     def _open_next_file(self):
#         self.cur_file_idx += 1
#         if self.cur_file_idx >= len(self.s5_files):
#             return False
#         logger.info(f"Switching to file {self.cur_file_idx} of {len(self.s5_files)} with name {self.s5_files[self.cur_file_idx]}")
#         self.cur_read_gen = get_slow5_reads_gen(str(self.s5_files[self.cur_file_idx]))
#         return True
    
#     def _fill_queue(self):
#         logger.info("Started queue filler thread")
#         while self._open_next_file():
#             for read in self.cur_read_gen:
#                 # logger.debug(f"Putting read {read[0]} into queue")
#                 self.queue.put(read)
#         self.queue.put(None) # sentinel
        
#     def get_new_read(self) -> Tuple[str, np.ndarray]:
#         # overwriting get_new_read directly rather than _get_new_read because queue already threadsafe
#         # res = self.queue.get()
#         # try getting an element immediately, otherwise print a warning
#         try:
#             res = self.queue.get_nowait()
#         except Empty:
#             logger.warning("Slow5 read queue empty, waiting until available")
#             res = self.queue.get()
#         if res is None:
#             self.definitely_empty = True
#             raise NoReadLeftException()
#         return res
    
#     # def stop_reading(self):
#     #     self.queue.put(None)
#     #     self._queue_filler_thread.join(), raise_if_error
    
#     # def reads_gen(self):
#     #     while True:
#     #         res = self.queue.get()
#     #         if res is None:
#     #             break
#     #         yield res
    
# s5_dir = "/home/mmordig/rawhash_project/rawhash2/test/data/d2_ecoli_r94/slow5_files"
# reader = Slow5ReadPool(s5_dir, 2)

# import time
# import tqdm

# num_signals = 0
# num_reads = 0
# # read_id, read_signal = reader.get_new_read()
# start_time = time.time()
# for (read_id, read_signal) in tqdm.tqdm(reader.reads_gen()):
#     num_reads += 1
#     num_signals += len(read_signal)
#     if num_reads > 10000:
#         break

# elapsed_time = time.time() - start_time
# num_signals / elapsed_time
# n_channels = 512
# min_throughput = n_channels * 4000
# num_signals / elapsed_time / min_throughput


# end_time_per_channel = np.array([time.time()] * n_channels)
# for (read_id, read_signal) in tqdm.tqdm(reader.reads_gen()):
#     min_i = np.argmin(end_time_per_channel)
#     time_sleep = end_time_per_channel[min_i] - time.time()
#     if time_sleep >= 0:
#         time.sleep(time_sleep)
#     # else:
#     #     print(f"Warning: missed deadline {time_sleep}")
#     end_time_per_channel[min_i] = time.time() + len(read_signal) / 4000
#     # xxx = read_signal.sum()