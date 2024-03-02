import pytest

from simreaduntil.shared_utils.nanosim_parsing import *

def test_parsing():
    
    assert NanoSimId.from_str("chr11-NC-000011_76599_perfect_proc0:0_F_0_9967_0").read_type == "perfect"
    
    # metagenomic read id
    nanosim_id_str = "Human-chr11-NC-000011_76599_aligned_proc0:0_F_0_9967_0"
    nanosim_id = NanoSimId.from_str(nanosim_id_str)
    assert nanosim_id.chrom == "Human-chr11-NC-000011"
    assert nanosim_id.ref_pos == 76599
    assert nanosim_id.read_nb == "proc0:0"
    assert nanosim_id.direction == "F"
    assert nanosim_id.head_len == 0
    assert nanosim_id.ref_len == 9967
    assert nanosim_id.tail_len == 0
    assert nanosim_id.read_type == "aligned"
    assert str(nanosim_id) == nanosim_id_str
    
    assert NanoSimId.is_valid(nanosim_id_str)
    assert not NanoSimId.is_valid("read1_765")
    
    with pytest.raises(AssertionError):
        # length larger than original
        nanosim_id.change_ref_len(10000)
        
    assert str(nanosim_id.change_ref_len(9967)) == nanosim_id_str # same length, so no m added
    
    assert str(nanosim_id.change_ref_len(1001)) == "Human-chr11-NC-000011_76599_aligned_proc0:0m_F_0_1001_0" # adds m to read_nb
    assert str(nanosim_id.change_ref_len(1001)) == "Human-chr11-NC-000011_76599_aligned_proc0:0m_F_0_1001_0" # same length, so no m added
    assert nanosim_id.ref_len == 1001 # in-place
    # reverse strand (R): 76599 + 9967 - 1001 = 85565
    assert str(NanoSimId.from_str("Human-chr11-NC-000011_76599_aligned_proc0:0_R_0_9967_0").change_ref_len(1001)) == "Human-chr11-NC-000011_85565_aligned_proc0:0m_R_0_1001_0" # adds m to read_nb
        
def test_parsing_unaligned():
    nanosim_id = NanoSimId.from_str("genome1-chr-6_236227_unaligned_proc5:16_R_0_16119_0")
    assert nanosim_id.chrom == "genome1-chr-6"
    assert nanosim_id.ref_pos == 236227
    assert nanosim_id.read_type == "unaligned"
    assert nanosim_id.read_nb == "proc5:16"
    assert nanosim_id.direction == "R"
    assert nanosim_id.head_len == 0
    assert nanosim_id.ref_len == 16119
    assert nanosim_id.tail_len == 0
    
def test_normalize_seq_name():
    assert normalize_seq_name("chr1 extra_info more-info hello") == "chr1-extra-info-more-info-hello"
    assert normalize_seq_name("chr1.aa extra_info more-info hello") == "chr1"

def test_case_convert_dna():
    assert case_convert_dna("AACCNNTGNAA", random_state=np.random.default_rng(2)) == 'AACCGTTGAAA'