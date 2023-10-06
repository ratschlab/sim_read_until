"""
Channel element within a channel, e.g., a read, a gap between reads, an unblock delay
"""

#%%
import enum
from numbers import Number
from functools import cached_property
from typing import Any, List, Optional, Tuple, Union
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId
from simreaduntil.shared_utils.numerics import really_less

#%%

class ChannelElement:
    """
    Base class for element/interval in channel [t_start, t_duration]
    
    Changing t_start or t_duration is possible.
    """
    def __init__(self, t_start: Number, t_duration: Number) -> None:
        self.t_start = t_start
        assert t_duration >= 0
        self.t_duration = t_duration
        
    @property
    def t_end(self):
        return self.t_start + self.t_duration
    
    @t_end.setter
    def t_end(self, t_end):
        assert t_end >= self.t_start
        self.t_duration = t_end - self.t_start
        
    def has_finished_by(self, t) -> bool:
        """Whether element has finished by time t"""
        return t >= self.t_end
    
    def elapsed_time(self, t: Number) -> Number:
        """Elapsed time of element up to time t, capped at element duration"""
        assert t >= self.t_start
        return min(t - self.t_start, self.t_duration)
    
    def __eq__(self, other) -> bool:
        """Equality, does not check elem_before, only time interval"""
        return self.t_start == other.t_start and self.t_duration == other.t_duration

def element_to_short_str(elem):
    """short string representation of element"""
    if isinstance(elem, ChunkedRead):
        return f"ChunkedRead({elem.full_read_id})"
    else:
        return type(elem).__name__
                       
class ShortGap(ChannelElement):
    """
    Short gap between reads in channel
    """
    def __init__(self, t_start, t_duration, elem_before: Optional["ChannelElement"]=None):
        super().__init__(t_start, t_duration)
        self.elem_before = elem_before
        
    def __repr__(self):
        res = f"ShortGap [{self.t_start}, {self.t_end}]"
        if self.elem_before is not None: 
            res += f" after element '{element_to_short_str(self.elem_before)}'"
        return res
    
class MuxScan(ChannelElement):
    """
    Mux scan in channel, so the channel is effectively unusable
    
    They are important to simulate ultra-long sequencing, where a long read may be broken due to a mux scan.
    
    Args:
        elem_to_restore: element to restore after mux scan
    """
    def __init__(self, t_start, t_duration, elem_to_restore: ChannelElement):
        super().__init__(t_start, t_duration)
        self.elem_to_restore = elem_to_restore
        
    def __repr__(self):
        res = f"MuxScan [{self.t_start}, {self.t_end}]"
        if self.elem_to_restore is not None: 
            res += f" with restorable element '{element_to_short_str(self.elem_to_restore)}'"
        return res
    
class UnblockDelay(ChannelElement):
    """
    Gap due to unblocking of read
    """
    def __init__(self, t_start, t_duration, read_before: Optional["ChannelElement"]=None):
        super().__init__(t_start, t_duration)
        self.read_before = read_before
        
    def __repr__(self):
        res = f"UnblockDelay [{self.t_start}, {self.t_end}]"
        if self.read_before is not None: 
            res += f" after read '{element_to_short_str(self.read_before)}'"
        return res
    
class LongGap(ChannelElement):
    """
    Long gap between reads due to temporary pore blockage (inactive period)
    """
    def __init__(self, t_start, t_duration, elem_before: Optional["ChannelElement"]=None):
        super().__init__(t_start, t_duration)
        self.elem_before = elem_before
        
    def split(self, t_split):
        """
        Split channel block into two at t_split, i.e. set self to [t_start, t_split] and return a new [t_split, t_end]
        
        Args:
            t_split: time to split at
            
        Returns:
            Second element
        """
        assert self.t_start <= t_split <= self.t_end
        orig_t_end = self.t_end
        self.t_duration = t_split - self.t_start
        return LongGap(t_split, orig_t_end - t_split, self.elem_before)
        
    def __repr__(self):
        res = f"PoreBlockage [{self.t_start}, {self.t_end}]"
        if self.elem_before is not None: 
            res += f" after element '{element_to_short_str(self.elem_before)}'"
        return res
    
class ChannelBroken(ChannelElement):
    """
    Infinite gap that indicates that the channel is completely broken and not just temporarily inactive
    """
    def __init__(self, t_start):
        super().__init__(t_start, np.inf)
        
    def __repr__(self):
        return f"ChannelBroken(since t={self.t_start})"

class NoReadLeftGap(ChannelElement):
    """
    Infinite gap because no more reads are left
    
    We use this class to avoid that the channel ever stops (i.e. should not throw an exception in a thread)
    """
    def __init__(self, t_start):
        super().__init__(t_start, np.inf)
        
    def __repr__(self):
        return f"NoReadLeftGap(since t={self.t_start})"

def estimate_ref_len(orig_ref_len, orig_seq_len, new_seq_len):
    """
    Estimate reference length of new sequence (cut short due to selective sequencing), not exact due to indels
    
    Args:
        orig_ref_len: original reference length
        orig_seq_len: original sequence length
        new_seq_len: new sequence length
        
    Returns:
        Estimated reference length of new sequence
    """
    return round(new_seq_len / orig_seq_len * orig_ref_len)

# StrEnum does not exist yet in Python3.8, see PythonDoc for IntEnum for this recipe
# allows printing as "field" instead of "class.field", where class is a class derived from enum
class ReadEndReason(str, enum.Enum):
    """Reason why a read ended"""
    
    # read ended prematurely
    SIM_STOPPED = "sim_stopped_unblocked" # simulation was stopped while read was still in-progress
    UNBLOCKED = "user_unblocked" # rejected read by ReadUntil
    MUX_SCAN_STARTED = "mux_scan_unblocked" # read was unblocked because a mux was started
    
    READ_ENDED_NORMALLY = "read_ended_normally" # read ended normally
    
end_reason_to_ont_map = {
    ReadEndReason.SIM_STOPPED.value: "unblock_mux_change",
    ReadEndReason.UNBLOCKED.value: "data_service_unblock_mux_change",
    ReadEndReason.MUX_SCAN_STARTED.value: "unblock_mux_change",
    ReadEndReason.READ_ENDED_NORMALLY.value: "signal_positive",
}
"""map simulator read end reasons to real ONT device read end reasons"""

class ReadTags(str, enum.Enum):
    """
    Tags to attach to a read, multiple are possible!
    """
    RU_NEVER_REQUESTED = "never_requested" # never requested by ReadUntil
    RU_STOPPED_RECEIVING = "stopped_receiving" # read was set to stop_receiving
    
class ChunkedRead(ChannelElement):
    """
    Read divided into chunks that can progressively be read from
    
    Chunks can be received from the read, it can be ended prematurely with .finish(), the sequence record or sequence summary record can be retrieved.
    If the read has a NanoSim read id, its id is modified to reflect the estimated reference length if it is ended prematurely.
    
    A read with n basepairs starts at t_start + t_delay, goes to time n*dt and the ith basepair (i>=1) is read after time i*dt, where dt=1/bp_per_second.
    Basepairs are returned in chunks of size chunk_size.
    
    t_start, t_delay, t_end, t_duration should not be modified once this class was instantiated. A read can be terminated early by calling .finish_now().
    
    Args:
        read_id: id of read
        seq: read sequence
        t_start: time at which the read starts
        read_speed: speed at which the read is read, in basepairs per second, defaults to SIM_PARAMS.bp_per_second
        chunk_size: size of chunks that .get_new_chunks() returns, defaults to SIM_PARAMS.chunk_size
        t_delay: delay before read starts (template_start - read_start), 0 basepairs are read during this delay, end time is shifted accordingly
    """
    def __init__(self, read_id: str, seq: str, t_start: Number, read_speed: Number=None, chunk_size: Number=None, t_delay:float = 0):
        # copy params in case they change while the read is in-progress
        assert read_speed > 0
        assert chunk_size > 0
        self._read_speed = read_speed
        self._chunk_size = chunk_size
        super().__init__(t_start,  len(seq) / self._read_speed + t_delay)
        
        self.full_read_id = read_id
        assert len(seq) > 0
        self._full_seq = seq
        assert t_delay >= 0
        self._t_delay = t_delay
        
        if NanoSimId.is_valid(read_id):
            # whether the read id is from NanoSim -> read id will be changed when read is terminated early
            self._nanosim_read_id = NanoSimId.from_str(read_id)
            self._ref_len = NanoSimId.from_str(self.full_read_id).ref_len
        else:
            self._nanosim_read_id = None
            # ref length not known, default to sequence length
            self._ref_len = len(self._full_seq)
        
        self.stopped_receiving = False # whether to receive chunks from the read
        self._next_chunk_idx = 0 # start idx of next chunks to return
        self.end_reason = None # action used to finish read
        
    def __repr__(self):
        return f"ChunkedRead '{self.full_read_id}': '{self._full_seq}' between [{self.t_start}, {self.t_end}], num chunks returned: {self._next_chunk_idx}, end_reason: {self.end_reason}"
    
    def __eq__(self, other) -> bool:
        assert isinstance(other, ChunkedRead)
        return self.t_start == other.t_start and self.t_duration == other.t_duration \
            and self._t_delay == other._t_delay \
            and self.full_read_id == other.full_read_id and self._full_seq == other._full_seq and self._read_speed == other._read_speed \
            and self._chunk_size == other._chunk_size and self.stopped_receiving == other.stopped_receiving \
            and self._next_chunk_idx == other._next_chunk_idx and self.end_reason == other.end_reason
    
    @property
    def _nb_chunks(self):
        """Number of chunks the read is divided into"""
        return (len(self._full_seq) + self._chunk_size - 1) // self._chunk_size # round up
    
    def full_duration(self) -> Number:
        return self._t_delay + len(self._full_seq) / self._read_speed
    
    @property
    def _has_received_chunks(self):
        """Whether at least one non-empty chunk was returned (via get_new_chunks)"""
        return self._next_chunk_idx > 0
    
    def all_chunks_consumed(self) -> bool:
        """Whether all chunks have been consumed, i.e. read has finished"""
        return self._next_chunk_idx >= self._nb_chunks
    
    @cached_property # lazy
    def _chunk_end_positions(self):
        """
        End positions of the chunks (cumulative chunk lengths), exclusive
        
        cached because they may not be needed if the simulation passes over the read (because time forwarded a lot).
        """
        cum_lens = [(i+1)*self._chunk_size for i in range(self._nb_chunks)]
        cum_lens[-1] = len(self._full_seq)
        return cum_lens
    
    def _get_chunks(self, fr: int, to: int):
        """Get concatenated chunks [from, to)"""
        return self._full_seq[fr * self._chunk_size:to * self._chunk_size]
    
    def _estimate_ref_len(self, nb_bps_read):
        """
        Estimate reference length given number of basepairs read, not exact due to indels
        
        Requires correct estimation of ref length of full read (which is available for NanoSim eads)
        """
        # round rather than round down (with int())
        assert nb_bps_read <= len(self._full_seq)
        return estimate_ref_len(orig_ref_len=self._ref_len, orig_seq_len=len(self._full_seq), new_seq_len=nb_bps_read)
    
    def nb_basepairs(self, t: Number):
        """
        Number of basepairs of read up to time t, first basepair emitted at time t_start
        
        If .has_finished_by(t) returns True and the full read was read or .finish() 
        not yet called, it is guaranteed to return the full length of the 
        read (independent of floating point errors).
        """
        real_start = self.t_start + self._t_delay
        if t < real_start:
            return 0
        if (self.has_finished_by(t) and self.end_reason in [None, ReadEndReason.READ_ENDED_NORMALLY]):
            # special case due to floating point problem with addition
            return len(self._full_seq)
        return min(
            len(self._full_seq), 
            int((min(t, self.t_end) - real_start) * self._read_speed) # round down
        )
        
    def nb_basepairs_full(self):
        """number of basepairs of full read"""
        return len(self._full_seq)
        
    def nb_basepairs_returned(self):
        """Number of basepairs not yet returned (if stopped receiving or if not getting chunks)"""
        return min(len(self._full_seq), self._next_chunk_idx * self._chunk_size)
        
    def finish(self, t=None, end_reason: Optional[ReadEndReason]=None) -> SeqIO.SeqRecord:
        """
        Finish read by time t, updating t_end
        
        This function can only be called once.
        
        Arguments:
            t: time when read ends, read is ended prematurely if less than full read end; full read if None or t exceeds full read end
            end_reason: action that caused read to be written to file
        Returns:
            SeqRecord of read that can be written to fasta file
        """
        assert self.end_reason is None, f"already ended read with action {self.end_reason}"
        
        if t is not None:
            if not self.has_finished_by(t):
                assert end_reason in [ReadEndReason.UNBLOCKED, ReadEndReason.MUX_SCAN_STARTED, ReadEndReason.SIM_STOPPED], "end reason must be set"
            
            nb_bps_returned = self.nb_basepairs_returned()
            assert t >= self.t_start + self._t_delay * (nb_bps_returned > 0) + nb_bps_returned / self._read_speed, "cannot finish earlier than last returned chunk"
            
            self.t_end = min(self.t_end, t)
        # t_end contains end time of read now
        
        if end_reason is None:
            end_reason = ReadEndReason.READ_ENDED_NORMALLY
        self.end_reason = end_reason
        
        return self.get_seq_record()
        
    def get_seq_record(self):
        """
        Get sequence record to write to FASTA file. Must call .finish() first.
        
        It contains extra information about read in the SeqRecord description field such as full read id, length.
        For NanoSim reads, the read id modified to reflect the actual read length if the read was terminated early.
        
        See ReadDescriptionParser for parsing the description field.
        """
        assert self.end_reason is not None, f"read was not ended yet"
        
        adapted_read_id = self.full_read_id
        if self.end_reason == ReadEndReason.READ_ENDED_NORMALLY:
            seq = self._full_seq
        else:
            seq = self._full_seq[:self.nb_basepairs(self.t_end)]
            if self._nanosim_read_id is not None and (len(seq) < len(self._full_seq)):
                # adapt reference length, as read was stopped early
                actual_ref_len = self._estimate_ref_len(self.nb_basepairs(self.t_end))
                self._nanosim_read_id.change_ref_len(actual_ref_len) # if this method is called again, the read id will not change again because the ref len is the same
                adapted_read_id = str(self._nanosim_read_id)
        
        # attach extra read tags
        read_tags = []
        if self.stopped_receiving:
            read_tags.append(ReadTags.RU_STOPPED_RECEIVING)
        if not self._has_received_chunks:
            read_tags.append(ReadTags.RU_NEVER_REQUESTED)
            
        # append full sequence length (in case read was unblocked)
        description = f"full_seqlen={self.nb_basepairs_full()} t_start={self.t_start} t_end={self.t_end} t_delay={self._t_delay} ended={self.end_reason} tags={','.join(read_tags)} full_read_id={self.full_read_id}"
        return SeqIO.SeqRecord(Seq(seq), id=adapted_read_id, description=description)
    
    
    SEQ_SUMMARY_HEADER = ["read_id", "channel", "mux", "start_time", "duration", "passes_filtering", "template_start", "template_duration", "sequence_length_template", "end_reason"]
    def get_seq_summary_record(self) -> Optional[List[str]]:
        """
        Get list of entries for sequence summary file, see SEQ_SUMMARY_HEADER for field names
        
        Returns:
            list of string entries, or None if read has no basepairs (e.g. due to delay)
        """
        # read_id, channel, mux, start_time, duration, passes_filtering, template_start, template_duration, sequence_length_template, end_reason
        mux = 1
        passes_filtering = "True"
        template_duration = self.t_duration - self._t_delay
        nb_bps_read = self.nb_basepairs(self.t_end)
        if nb_bps_read <= 0:
            return None
        return [
            self.full_read_id, self.channel, mux, self.t_start, self.t_duration, passes_filtering, 
            self.t_start + self._t_delay, template_duration, 
            nb_bps_read, self.end_reason
        ]
        
    def stop_receiving(self, value=True) -> bool:
        """
        Set the read to stop receiving, so no new chunks will be received from this read (with get_new_chunk), but it keeps being read
        
        Arguments:
            value: whether to stop receiving
        Returns:
            True iff the read's stop_receiving state changed
        """
        assert self.end_reason is None, f"already ended read with action {self.end_reason}"
        
        if self.stopped_receiving == value:
            return False
        self.stopped_receiving = value
        return True
    
    def get_new_chunks(self, t: Number) -> Tuple[str, str, Optional[int]]:
        """
        Get new read chunks up to time <= t, only new data since last call to this function
        
        Implicitly forwards to time t (choosing chunk index of chunk containing t)
        
        Returns:
            (all chunks concatenated, read_id, estimated_ref_len_so_far).
            
            Empty chunks "" if stop_receiving is True
            estimated_ref_len_so_far is the estimated number of basepairs covered by chunks returned so far
        """
        assert self.end_reason is None, f"already ended read with action {self.end_reason}"
        
        if self.stopped_receiving:
            return "", self.full_read_id, self._estimate_ref_len(self.nb_basepairs_returned())
        
        # choose chunk index of chunk containing t
        next_chunk_idx = np.searchsorted(self._chunk_end_positions, v=self.nb_basepairs(t), side='right') # index such that a[i-1] <= v < a[i]
        old_next_chunk_idx = self._next_chunk_idx
        self._next_chunk_idx = next_chunk_idx
        estimated_ref_len_so_far = self._estimate_ref_len(self.nb_basepairs_returned()) # takes into account _next_chunk_idx
        return self._get_chunks(old_next_chunk_idx, self._next_chunk_idx), self.full_read_id, estimated_ref_len_so_far


class ReadDescriptionParser:
    """
    Parse the read description field written by the ChunkedRead
    
    The description must have the right fields (as the simulator ensures), otherwise this function will fail.
    
    Args:
        description: the description field of a SeqRecord; pay attention to SeqIO.SeqRecord.description, which contains the id 
            in the description as well when parsed from a file, so do record.description.split(" ", maxsplit=1)[1] to get the description
    """
    def __init__(self, description: str):
        self._fields = {} # set it, otherwise overriden getattribute fails
        self._parse_description(description)
        
    def _parse_description(self, description: str):
        elements = description.split(" ")
        tags = dict((elem.split("=") for elem in elements))
        tags["full_seqlen"] = int(tags["full_seqlen"])
        tags["t_start"] = float(tags["t_start"])
        tags["t_end"] = float(tags["t_end"])
        tags["t_delay"] = float(tags["t_delay"])
        tags["tags"] = tags["tags"].split(",")
        self._fields = tags
        
        self._description = description
        
    def __repr__(self):
        return f"ReadDescriptionParser({self._description}, {self._fields})"
    
    def __getattribute__(self, name: str):
        """
        Make tags (saved in self._fields) available as attributes, e.g. self.full_seqlen
        """
        if name == "_fields" or name not in self._fields.keys():
            return object.__getattribute__(self, name)
        return self._fields[name]
