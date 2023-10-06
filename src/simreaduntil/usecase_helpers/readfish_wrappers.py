"""
Wrappers to make the ONTSimulator work with ReadFish
"""

from contextlib import contextmanager
import logging
import time
from typing import Iterable, Tuple

import numpy as np
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId
from simreaduntil.shared_utils.timing import cur_ns_time

from simreaduntil.simulator.simulator import ReadUntilDevice
from ru.utils import Severity

logger = setup_logger_simple(__name__)
"""module logger"""

class ReadWrapper:
    """
    Wrapper around read for .get_read_chunks().
    
    read_nb and read_id are the same thing.
    """
    def __init__(self, read_id, seq: str, quality_seq: str) -> None:
        self.id = read_id
        assert isinstance(seq, str)
        self.seq = seq
        self.quality_seq = quality_seq
        
    @property
    def number(self):
        return self.id
    
    def __repr__(self) -> str:
        return f"ReadWrapper({self.id}, {self.seq}, {self.quality_seq})"
    
class ReadUntilClientWrapper:
    """
    Wraps a ReadUntilDevice and provides a read_until interface to substitute ReadFish's ru.read_until_client
    
    The read chunks are already basecalled, not raw signal data.
    read_nb and read_id are the same thing.
    
    Args:
        client: ReadUntilDevice
    """
    def __init__(self, client: ReadUntilDevice) -> None:
        self._client = client
        self.connection = None
        
    def __repr__(self) -> str:
        return f"ReadUntilClientWrapper({self._client.__repr__()})"
    
    @property
    def is_running(self):
        return self._client.is_running
    
    is_phase_sequencing = is_running
    signal_dtype = None
    
    @property
    def mk_run_dir(self):
        return self._client.mk_run_dir
    
    def unblock_read_batch(self, unblock_batch_action_list, duration):
        self._client.unblock_read_batch([(channel, read_id) for (channel, read_nb, read_id) in unblock_batch_action_list], unblock_duration=duration)
    
    def stop_receiving_batch(self, stop_receiving_action_list):
        self._client.stop_receiving_batch(stop_receiving_action_list)
    
    def get_read_chunks(self, batch_size=1, last=True):
        for (channel, read_id, seq, quality_seq, estimated_ref_len_so_far) in self._client.get_basecalled_read_chunks(batch_size=batch_size):
            yield (channel, ReadWrapper(read_id=read_id, seq=seq, quality_seq=quality_seq))


BASECALLING_TIME_PER_BP = (4000/450) / 2.6e7
"""
Time to basecall one basepair

.. python::
    # 30Gbps in 72h for MinION, see https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview
    # this accounts for channel deterioration, so the bps/s is lower than at peak
    bps_per_second_lower = 30 * 1e9 / (72 * 24 * 3600)
    print(f"Average basecalling speed: {bps_per_second_lower:e} bp/s")

    Another estimate is to assume 80% of the pores are reading (more at the beginning, much less in the end)
    1 / (450 * 512 * 0.8) = 5.4e-6 s/bp
    This assumes that the basecaller has just the right speed, though it may actually be faster in practice, as considered below.

    # see https://hackmd.io/@Miles/HJUnkIeOK for Tesla V100 (similar to GV100 GPU in GridION)
    # 2.6e7 samples/s where sample = 1 current signal, 4kHz, 450bp/s, so 4000/450 = 8.9 samples/bp
    bps_per_second = 2.6e7 / (4000/450) = 3.4e-7 s/bp
    time_per_bp = 1 / bps_per_second
    print(f"Basecalling speed: {bps_per_second:e} bp/s, time per basepair {time_per_bp:e} s/bp")

"""

class DummyBasecaller:
    """
    Dummy basecaller not doing anything, just returning already basecalled sequences
    
    Imitating the GuppyCaller class from ReadFish
    Since basecalling is not done, the basecaller sleeps for some time.
    
    Args:
        time_per_bp: time to basecall one basepair
    """
    def __init__(self, time_per_bp=None):
        if time_per_bp is None:
            time_per_bp = BASECALLING_TIME_PER_BP
        assert time_per_bp >= 0, f"Basecalling speed {time_per_bp} is negative"
        self.time_per_bp = time_per_bp
    
    def basecall_minknow(self, reads: Iterable[Tuple[int, ReadWrapper]], signal_dtype: str, decided_reads: Iterable[Tuple[int, int]]):
        """
        Basecall data from minknow

        Args:
            reads: List or generator of tuples containing (channel, ReadWrapper)
            signal_dtype: ignored
            decided_reads: ignored
        
        Yields:
            (channel, read_number), read_id, sequence, sequence_length, quality
        """
        
        time_start = time.perf_counter_ns() # in nanoseconds, only offsets are correct
        total_wait_time = 0
        nb_called_bps = 0
        for (channel, read_info) in reads:
            # if self.time_per_bp > 0:
            #     time.sleep(len(read_info.seq) * self.time_per_bp)
            # to imitate the guppy basecaller which runs in parallel, we do not delay each time something is requested, but rather since the function was called
            nb_called_bps += len(read_info.seq)
            if self.time_per_bp > 0:
                wait_time = nb_called_bps * self.time_per_bp - (time.perf_counter_ns() - time_start)/1_000_000_000
                if wait_time > 0:
                    time.sleep(wait_time)
                    total_wait_time += wait_time
            yield (channel, read_info.number), read_info.id, read_info.seq, len(read_info.seq), read_info.quality_seq
            
        if nb_called_bps > 0:
            logger.info(f"Total basecalling extra wait time: {total_wait_time:.2e}s for {nb_called_bps} basepairs")

device_logger = logging.getLogger("device_logger")
"""logger to use for messages to the ONT device"""

def send_message_to_logger(rpc_connection, message, severity: Severity):
    """
    Send a message to a logger rather than to the real device
    
    You have to replace the read_until.send_message function with this one
    """
    log_method = None
    if severity == Severity.INFO:
        log_method = device_logger.info
    elif severity == Severity.WARN:
        log_method = device_logger.warning
    else:
        assert severity == Severity.ERROR
        # logger.error does not show traceback by default
        log_method = device_logger.exception
        
    log_method(f"DEVICE MESSAGE: {message} (severity {severity})")
    
def get_flowcell_array_replacement(flowcell_size):
    """
    Replacement for ru.get_flowcell_array to work with any flowcell size
    """
    # return np.array
    return np.arange(1, flowcell_size + 1).reshape(1, -1)

class NanoSimMapper:
    """
    Align NanoSim read by parsing its id
    
    This class avoids minimap2, so the mapping delay can be controlled when running in accelerated mode
    """
    class Alignment:
        def __init__(self, query_name, query_len, query_start, query_end, target_strand, target_name, target_len, target_start, target_end, num_matches, alignment_block_length, mapping_quality):
            self.query_name = query_name
            self.query_len = query_len
            self.query_start = query_start
            self.query_end = query_end
            self.strand = target_strand
            self.ctg = target_name
            self.target_len = target_len
            self.r_st = target_start
            self.target_end = target_end
            self.num_matches = num_matches
            self.alignment_block_length = alignment_block_length
            self.mapping_quality = mapping_quality
    
        def __repr__(self):
            return f"{self.query_name}\t{self.query_len}\t{self.query_start}\t{self.query_end}\t{self.strand}\t{self.ctg}\t{self.target_len}\t{self.r_st}\t{self.target_end}\t{self.num_matches}\t{self.alignment_block_length}\t{self.mapping_quality}"
        
    def __init__(self, index):
        self.initialised = True
        self.index = index # not used
        
    @staticmethod
    def _map_seq(read_id, seq_len):
        parsed = NanoSimId.from_str(read_id)
        
        return NanoSimMapper.Alignment(
            query_name=read_id, query_len=seq_len, query_start=0, query_end=seq_len, 
            target_strand=1 if parsed.direction == "F" else -1, target_name=parsed.chrom, target_len="*", target_start=parsed.ref_pos, target_end=parsed.ref_len, 
            num_matches=seq_len, alignment_block_length=seq_len, mapping_quality=255
        )
    
    def map_reads_2(self, calls):
        """Align reads against a reference

        Args:
            calls: iterable of called reads from PerpetualCaller.basecall_minknow, iterable [tuple,  str, str, int, str]

        Yields
            tuple ((channel, read_number), read_id, sequence, sequence_length, mapping_results)
        """
        for read_info, read_id, seq, seq_len, quality in calls:
            assert len(seq) == seq_len
            yield read_info, read_id, seq_len, [self._map_seq(read_id, seq_len)]
            
@contextmanager
def replace_ru_mapper(replace):
    """
    Contextmanager to replace ru mapper with NanoSim mapper
    
    Args:
        replace: whether to replace it
        
    Yields:
        mapper class to use
    """
    import ru.ru_gen
    original_mapper = ru.ru_gen.CustomMapper
    if replace:
        logger.info("Using fake minimap2 mapper")
        
        # cannot replace ru.basecall.Mapper since ru.ru_gen imports it directly
        ru.ru_gen.CustomMapper = NanoSimMapper
        
        yield NanoSimMapper
        
        # restore
        ru.ru_gen.CustomMapper = original_mapper
    else:
        yield original_mapper
        