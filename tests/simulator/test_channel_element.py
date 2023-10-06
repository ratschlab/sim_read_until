import numpy as np
import pytest
from pytest import approx

from Bio import SeqIO
from Bio.Seq import Seq

from simreaduntil.simulator.channel_element import ShortGap, ChunkedRead as _ChunkedRead, NoReadLeftGap, LongGap, ReadDescriptionParser, ReadEndReason, end_reason_to_ont_map

ChunkedRead = lambda *args, **kwargs: _ChunkedRead(*args, **kwargs, read_speed=10, chunk_size=4)
    
eps = 1e-8 # small delay (to avoid issues with rounding errors when geting chunks up to time <= t)
    
def test_readended_map():
    assert end_reason_to_ont_map[ReadEndReason.UNBLOCKED.value] == "data_service_unblock_mux_change"
    
def test_nb_basepairs():
    # test basic functions
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1)
    assert chunked_read._nanosim_read_id is None
    assert chunked_read._nb_chunks == 5
    assert chunked_read._chunk_end_positions == [4, 8, 12, 16, 18]
    assert chunked_read._get_chunks(2, 4) == "33334444"
    assert chunked_read._get_chunks(4, 10) == "55"
    assert chunked_read.nb_basepairs(10.1+eps) == 0
    assert chunked_read.nb_basepairs(10.18+eps) == 0
    assert chunked_read.nb_basepairs(10.2+eps) == 1
    assert chunked_read.nb_basepairs(10.35+eps) == 2
    assert chunked_read.nb_basepairs(13) == len(chunked_read._full_seq)
    assert chunked_read.t_end == approx(10.1 + 1.8)
    assert chunked_read.full_duration() == approx(1.8)
    
    # has_finished
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1)
    assert not chunked_read.has_finished_by(10.3)
    assert chunked_read.has_finished_by(20.3)
    
    chunked_read.finish(10.9+eps, end_reason=ReadEndReason.UNBLOCKED)
    assert chunked_read.t_end == 10.9+eps
    assert chunked_read.nb_basepairs(10.9+eps) == 8
    assert chunked_read.nb_basepairs_full() == 18
    
    nanosim_id = "chr11_77_aligned_proc0:0_F_0_36_0"
    chunked_read = ChunkedRead(nanosim_id, "111122223333444455", 10.1)
    assert str(chunked_read._nanosim_read_id) == nanosim_id
    
    # check that all basepairs are returned once the read has finished
    # issue previously due to floating point
    chunked_read = ChunkedRead("read1", "111122223", 28.4)
    t_end = 29.299999999999997
    chunked_read.has_finished_by(t_end)
    assert chunked_read.nb_basepairs(t_end) == 9
    
def test_chunks():
    # check get_new_chunks
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1)
    assert not chunked_read.all_chunks_consumed()
    # 1 bp emitted every 0.1 seconds, add small tolerance if it is just on the edge
    assert chunked_read.get_new_chunks(8.9) == ("", "read1", 0)
    assert chunked_read.nb_basepairs_returned() == 0
    assert chunked_read.get_new_chunks(10.1+0.3+eps) == ("", "read1", 0)
    assert chunked_read.get_new_chunks(10.1+0.4+eps) == ("1111", "read1", 4)
    assert chunked_read.nb_basepairs_returned() == 4
    assert chunked_read.get_new_chunks(10.1+0.4+eps) == ("", "read1", 4), "no new chunks since last time"
    assert not chunked_read.all_chunks_consumed()
    assert chunked_read.get_new_chunks(10.1+0.4+eps) == ("", "read1", 4), "no new chunks since last time"
    assert chunked_read.get_new_chunks(10.1+0.4+1+eps) == ("22223333", "read1", 12)
    assert chunked_read.nb_basepairs_returned() == 12
    assert chunked_read.get_new_chunks(10.1+0.4+1.+0.7+eps) == ("444455", "read1", 18)
    assert chunked_read.get_new_chunks(130.2) == ("", "read1", 18)
    assert chunked_read.nb_basepairs_returned() == 18
    assert chunked_read.all_chunks_consumed()
    
    # test with NanoSim-like read, ref_len = 36, start position = 77
    # 36/18 = 2
    nanosim_read_id = "chr11_77_aligned_proc0:0_F_0_36_0"
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 10.1)
    assert chunked_read.get_new_chunks(10.1+0.9+eps) == ("11112222", nanosim_read_id, 16)
    assert chunked_read.get_new_chunks(10.1+3.9+eps) == ("3333444455", nanosim_read_id, 36)
    
    # test stop_receiving
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1)
    chunked_read.stop_receiving()
    assert chunked_read.get_new_chunks(10.1+0.9+eps) == ("", "read1", 0)
    
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1)
    assert chunked_read.get_new_chunks(10.1+0.9+eps) == ("11112222", "read1", 8)
    chunked_read.stop_receiving()
    assert chunked_read.get_new_chunks(10.1+1.4+eps) == ("", "read1", 8)
    
def test_chunks_with_delay():
    # extra delay before actual read starts
    chunked_read = ChunkedRead("read1", "111122223333444455", 10.1, t_delay=2.1)
    assert chunked_read._nb_chunks == 5
    assert chunked_read._chunk_end_positions == [4, 8, 12, 16, 18]
    assert chunked_read.nb_basepairs(10.4) == 0
    assert chunked_read.nb_basepairs(10.1 + 2.15) == 0
    assert chunked_read.nb_basepairs(10.1 + 2.25) == 1
    assert chunked_read.nb_basepairs(10.1 + 3.25) == 11
    assert chunked_read.t_end == approx(10.1 + 2.1 + 1.8)
    assert chunked_read.full_duration() == approx(1.8 + 2.1)
    
    assert chunked_read.get_new_chunks(8.9) == ("", "read1", 0)
    assert chunked_read.get_new_chunks(10.5) == ("", "read1", 0)
    assert chunked_read.get_new_chunks(10.1 + 2.3) == ("", "read1", 0)
    assert chunked_read.get_new_chunks(10.1+2.1+0.4+eps) == ("1111", "read1", 4)
    assert chunked_read.get_new_chunks(10.1+5+eps) == ("22223333444455", "read1", 18)
    
    # test finish
    nanosim_read_id = "chr11_77_aligned_proc0:0_F_0_36_0"
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2, t_delay=2.1)
    assert chunked_read.get_new_chunks(0.2+2.1+0.9+eps) == ("11112222", nanosim_read_id, 16) # 2 * 8 = 16
    with pytest.raises(AssertionError, match="finish earlier"):
        chunked_read.finish(0.2+2.1+0.5+eps, end_reason=ReadEndReason.UNBLOCKED)
    
    seq_record = chunked_read.finish(0.2+2.1+1.0+eps, end_reason=ReadEndReason.UNBLOCKED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_77_aligned_proc0:0m_F_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={0.2+2.1+1.0+eps} t_delay={2.1} ended=user_unblocked tags= full_read_id={nanosim_read_id}")) # 10*2 = 20
    
def check_equal_seq_records(seq_record1, seq_record2):
    # using assert here rather than return a boolean has the advantage that pytest gives good error messages (which inequality failed)
    assert seq_record1.id == seq_record2.id and seq_record1.seq == seq_record2.seq and seq_record1.description == seq_record2.description

def test_read_finish():
    # starting at position 77 and ref length 36, 36/18=2
    nanosim_read_id = "chr11_77_aligned_proc0:0_F_0_36_0"
    
    # get chunks, but no action
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish()
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("111122223333444455"), id=nanosim_read_id, description=f"full_seqlen=18 t_start=0.2 t_end={0.2+1.8} t_delay={0} ended=read_ended_normally tags= full_read_id={nanosim_read_id}"))
    assert chunked_read.end_reason == ReadEndReason.READ_ENDED_NORMALLY
    with pytest.raises(AssertionError, match="already ended"):
        # cannot finish again
        chunked_read.finish()
    
    # negative start time, read finished
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", -0.5)
    chunked_read.get_new_chunks(0.9)
    chunked_read.finish(-0.5+1.9, end_reason=ReadEndReason.READ_ENDED_NORMALLY)
    
    # get chunks, then stop receiving
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    chunked_read.stop_receiving()
    seq_record = chunked_read.finish()
    check_equal_seq_records(seq_record, chunked_read.get_seq_record())
    check_equal_seq_records(seq_record, chunked_read.get_seq_record()) # check idempotency
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("111122223333444455"), id=nanosim_read_id, description=f"full_seqlen=18 t_start=0.2 t_end={0.2+1.8} t_delay={0} ended=read_ended_normally tags=stopped_receiving full_read_id={nanosim_read_id}"))
    assert chunked_read.end_reason == ReadEndReason.READ_ENDED_NORMALLY
    
    # never requested chunks, but still stop receiving
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.stop_receiving()
    seq_record = chunked_read.finish()
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("111122223333444455"), id=nanosim_read_id, description=f"full_seqlen=18 t_start=0.2 t_end={0.2+1.8} t_delay={0} ended=read_ended_normally tags=stopped_receiving,never_requested full_read_id={nanosim_read_id}"))
    assert chunked_read.end_reason == ReadEndReason.READ_ENDED_NORMALLY
    
    # never requested chunks
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    seq_record = chunked_read.finish()
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("111122223333444455"), id=nanosim_read_id, description=f"full_seqlen=18 t_start=0.2 t_end={0.2+1.8} t_delay={0} ended=read_ended_normally tags=never_requested full_read_id={nanosim_read_id}"))
    assert chunked_read.end_reason == ReadEndReason.READ_ENDED_NORMALLY
    
    # stopped receiving, rejected afterwards
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    chunked_read.stop_receiving()
    chunked_read.get_new_chunks(1.1)
    seq_record = chunked_read.finish(1.2+eps, end_reason=ReadEndReason.UNBLOCKED)
    check_equal_seq_records(seq_record, chunked_read.get_seq_record())
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_77_aligned_proc0:0m_F_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={1.2+eps} t_delay={0} ended=user_unblocked tags=stopped_receiving full_read_id={nanosim_read_id}")) # 10*2 = 20
    assert chunked_read.end_reason == ReadEndReason.UNBLOCKED
    
    # sim stopped
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish(1.2+eps, end_reason=ReadEndReason.SIM_STOPPED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_77_aligned_proc0:0m_F_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={1.2+eps} t_delay={0} ended=sim_stopped_unblocked tags= full_read_id={nanosim_read_id}")) # 10*2 = 20
    assert chunked_read.end_reason == ReadEndReason.SIM_STOPPED
    
    # mux scan started
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish(1.2+eps, end_reason=ReadEndReason.MUX_SCAN_STARTED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_77_aligned_proc0:0m_F_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={1.2+eps} t_delay={0} ended=mux_scan_unblocked tags= full_read_id={nanosim_read_id}")) # 10*2 = 20
    assert chunked_read.end_reason == ReadEndReason.MUX_SCAN_STARTED
    
    # terminate early without end reason
    with pytest.raises(AssertionError, match="end reason"):
        chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
        chunked_read.get_new_chunks(0.9)
        chunked_read.finish(1.1) # need to indicate end reason
        
    # cannot finish earlier than last chunk received
    with pytest.raises(AssertionError, match="finish earlier"):
        chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
        chunked_read.get_new_chunks(0.9)
        chunked_read.finish(0.5, end_reason=ReadEndReason.UNBLOCKED)
        
    # can finish at 0.6 since last chunk returned then
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish(0.6+eps, end_reason=ReadEndReason.UNBLOCKED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111"), id="chr11_77_aligned_proc0:0m_F_0_8_0", description=f"full_seqlen=18 t_start=0.2 t_end={0.6+eps} t_delay={0} ended=user_unblocked tags= full_read_id={nanosim_read_id}")) # 4*2 = 8
    assert chunked_read.end_reason == ReadEndReason.UNBLOCKED
        
    # unblock read
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish(1.2+eps, end_reason=ReadEndReason.UNBLOCKED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_77_aligned_proc0:0m_F_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={1.2+eps} t_delay={0} ended=user_unblocked tags= full_read_id={nanosim_read_id}")) # 10*2 = 20
    assert chunked_read.end_reason == ReadEndReason.UNBLOCKED
    
    # unblock on reverse strand
    # starting at position 77 and ref length 36, 36/18=2
    nanosim_read_id = "chr11_77_aligned_proc0:0_R_0_36_0"
    chunked_read = ChunkedRead(nanosim_read_id, "111122223333444455", 0.2)
    chunked_read.get_new_chunks(0.9)
    seq_record = chunked_read.finish(1.2+eps, end_reason=ReadEndReason.UNBLOCKED)
    check_equal_seq_records(seq_record, SeqIO.SeqRecord(Seq("1111222233"), id="chr11_93_aligned_proc0:0m_R_0_20_0", description=f"full_seqlen=18 t_start=0.2 t_end={1.2+eps} t_delay={0} ended=user_unblocked tags= full_read_id={nanosim_read_id}")) # 77 + 36 - 10*2 = 93
    assert chunked_read.end_reason == ReadEndReason.UNBLOCKED
    
def test_no_readleft_gap():
    no_read_left = NoReadLeftGap(3)
    assert not no_read_left.has_finished_by(2)
    assert not no_read_left.has_finished_by(3+eps)
    assert not no_read_left.has_finished_by(5+eps)
    
def test_ChannelGap():
    elem_before = ShortGap(0, 1) # end time does not matter
    gap = LongGap(3, t_duration=7, elem_before=elem_before)
    next_gap = gap.split(8)
    assert gap.t_start == 3
    assert gap.t_end == 8
    assert gap.elem_before is elem_before
    assert next_gap.t_start == 8
    assert next_gap.t_end == 10
    assert next_gap.elem_before is elem_before
    
    with pytest.raises(AssertionError):
        LongGap(3, t_duration=7).split(12)
        
def test_ReadDescriptionParser():
    description = "full_seqlen=18 t_start=0.2 t_end=2 t_delay=0.1 ended=read_ended_normally tags=stopped_receiving,never_requested"
    parsed_desc = ReadDescriptionParser(description)
    assert parsed_desc._fields == {'full_seqlen': 18, 't_start': 0.2, 't_end': 2.0, 't_delay': 0.1, 'ended': 'read_ended_normally', 'tags': ['stopped_receiving', 'never_requested']}
    assert parsed_desc.full_seqlen == 18
    print(parsed_desc)

    # empty tags
    description = "full_seqlen=18 t_start=0.2 t_end=2 t_delay=0.1 ended=read_ended_normally tags= otherfield=othervalue"
    ReadDescriptionParser(description)