"""
Read pool that returns reads when requested, e.g., from a generator or a file
"""

import threading
from typing import Optional, Dict, Any, Generator
import numpy as np
import pysam

from simreaduntil.shared_utils.utils import force_eval_generator_function, get_some_value_from_dict, is_empty_file

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
    
class NoReadLeft(Exception):
    """
    When no read is left in the read pool
    """
    pass

class ReadPool:
    """
    Read pool from which reads can be obtained
    """
    def __init__(self):
        self.lock = threading.Lock()
        self.definitely_empty = False # whether the readpool is definitely empty (if False, we don't know)
        self.nb_reads_returned = 0
        
    def get_new_read(self, channel=None) -> str:
        """
        Get new read (thread-safe)
        
        Once the function returns NoReadLeft for a channel, it will always return NoReadLeft for this channel (also for channel=None).
        
        Args:
            channel: channel for which to get read
        """
        with self.lock:
            res = self._get_new_read(channel=channel)
            self.nb_reads_returned += 1
            return res
    
    def _get_new_read(self, channel=None) -> str:
        """
        Get new read (not thread-safe), overwrite in subclasses
        
        Args:
            channel: channel for which to get read
        """
        raise NotImplementedError()

class ReadPoolFromIterable(ReadPool):
    """
    Read pool that requests reads from generator
    """
    def __init__(self, reads_iterable):
        super().__init__()
        
        self.reads_iterable = reads_iterable
        
    def _get_new_read(self, channel=None) -> str:
        # note: generators are not thread-safe !!
        try:
            return next(self.reads_iterable)
        except StopIteration as e:
            self.definitely_empty = True
            raise NoReadLeft from e # exception chaining
        
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
        super().__init__()
        
        self.reads_iterable_per_channel = reads_iterable_per_channel
        
    def _get_new_read(self, channel=None) -> str:
        try:
            return next(self.reads_iterable_per_channel[channel])
        except StopIteration as e:
            raise NoReadLeft from e # exception chaining
        
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
    Keep track of reads_file
    """
    def __init__(self, reads_file, shuffle_rand_state: Optional[np.random.Generator]=None):
        super().__init__(reads_from_file_gen(reads_file, shuffle_rand_state=shuffle_rand_state))
        
        self.shuffled = shuffle_rand_state is not None
        self.reads_file = reads_file
        
    def __repr__(self):
        return f"ReadPool(file = {self.reads_file}, shuffled = {self.shuffled})"
    
    