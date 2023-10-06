import copy
from textwrap import dedent
from matplotlib.layout_engine import TightLayoutEngine
import numpy as np
import pytest
from pytest import approx
from simreaduntil.shared_utils.testing_utils import assert_dict_with_np_arrays_equal

from simreaduntil.seqsum_tools.coverage_tracker import NanoSimCoverageTracker, PafCoverageTracker, _add_count_to_blocked_array, _divide_round_up, plot_cov_per_position_lineplot, plot_covs_hist

# todo2: remove blocksize, CovTracker.empty_from_lens("nanosim", blocksize=10, *args): i.e. don't specify blocksize, but compute good value for it
def test_plot_cov_per_bp_lineplot():
    covs_per_bp = np.array([1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 3, 3, 0, 0, 5, 5, 5, 5, 1, 1, 1, 0, 0, 0, 0, 0])
    plot_cov_per_position_lineplot(covs_per_bp, n_points=100, figsize=(10, 5))
    plot_cov_per_position_lineplot(covs_per_bp, n_points=100, figsize=(10, 5), target_coverage=4)
    
def test_plot_covs_hist():
    covs_per_bp = np.array([1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 3, 3, 0, 0, 5, 5, 5, 5, 1, 1, 1, 0, 0, 0, 0, 0])
    plot_covs_hist(covs_per_bp)
    plot_covs_hist(covs_per_bp, target_coverage=4)
    
def test_divide_round_up():
    assert _divide_round_up(0, 1) == 0
    assert _divide_round_up(1, 1) == 1
    assert _divide_round_up(2, 1) == 2
    assert _divide_round_up(0, 3) == 0
    assert _divide_round_up(1, 3) == 1
    assert _divide_round_up(3, 3) == 1
    assert _divide_round_up(4, 3) == 2
    
def test_add_count_to_blocked_array():
    def add_count_to_empty(*args, **kwargs):
        arr = np.array([0, 0, 0, 0, 0], dtype=np.uint8)
        _add_count_to_blocked_array(arr, *args, **kwargs)
        return arr

    assert add_count_to_empty(0, 2, 3) == approx(np.array([2, 0, 0, 0, 0]))
    assert add_count_to_empty(0, 4, 3) == approx(np.array([3, 1, 0, 0, 0]))
    assert add_count_to_empty(2, 7, 3) == approx(np.array([1, 3, 1, 0, 0]))
    assert add_count_to_empty(4, 10, 3) == approx(np.array([0, 2, 3, 1, 0]))
    assert add_count_to_empty(4, 10, 3) == approx(np.array([0, 2, 3, 1, 0]))
    assert add_count_to_empty(4, 15, 3) == approx(np.array([0, 2, 3, 3, 3]))
    assert add_count_to_empty(14, 15, 3) == approx(np.array([0, 0, 0, 0, 1]))
    assert add_count_to_empty(4, 4, 3) == approx(np.array([0, 0, 0, 0, 0]))

    assert add_count_to_empty(0, 5, 1) == approx(np.array([1, 1, 1, 1, 1]))
    assert add_count_to_empty(2, 4, 1) == approx(np.array([0, 0, 1, 1, 0]))
    assert add_count_to_empty(2, 4, 1) == approx(np.array([0, 0, 1, 1, 0]))

    random_state = np.random.default_rng(2)
    blocksize = 7
    cov_arr = np.zeros((300 + blocksize-1) // blocksize, dtype=np.int32)
    expected_sum = 0
    for i in range(100): # < 256 -> fits into count
        start = random_state.integers(0, 300)
        length = random_state.integers(1, 300 - start + 1)
        cov_arr_prev = cov_arr.copy()
        _add_count_to_blocked_array(cov_arr, start, start + length, blocksize=blocksize)
        assert all(cov_arr >= cov_arr_prev)
        expected_sum += length
        if cov_arr.sum() != expected_sum:
            print(f"{i}: {start}, {start+length}")
        assert cov_arr.sum() == expected_sum
        
##### CovTracker without blocks
def test_covtracker_from_empty(tmp_path):
    ref_genome_file = tmp_path / "ref_genome.fa"
    with open(ref_genome_file, "w") as f:
        f.write(dedent("""\
            > chr1
            AAATTTTGGGTTTCCCCA
            > chr2
            CCTTAAATTAAAATTCCCCCCCGGGGG
            """))
        
    assert_dict_with_np_arrays_equal(NanoSimCoverageTracker.empty_from_ref_genomes(ref_genome_file, blocksize=1).coverage_per_chrom, {"chr1": np.zeros(18), "chr2": np.zeros(27)})
    
    assert_dict_with_np_arrays_equal(NanoSimCoverageTracker.empty_from_lens({"chr1": 18, "chr2": 27}, blocksize=1).coverage_per_chrom, {"chr1": np.zeros(18), "chr2": np.zeros(27)})

def test_covtracker_math():
    # test functions that compute something
    covs_per_bp_per_chrom = {
        "chr1": np.array([1, 1, 2, 3, 2]),
        "chr2": np.array([3, 3, 2, 1, 1, 2, 2, 2]),
        "chr3": np.array([1, 1, 3, 3])
    }
    cov_tracker = NanoSimCoverageTracker(copy.deepcopy(covs_per_bp_per_chrom), {chrom: len(covs) for chrom, covs in covs_per_bp_per_chrom.items()}, blocksize=1)
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr1"]) == approx(3/5)
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr2"]) == approx(6/8)
    assert cov_tracker.get_fraction_cov_atleast(3, ["chr1"]) == approx(1/5)
    assert cov_tracker.get_fraction_cov_atleast(3, ["chr3"]) == approx(2/4)
    assert cov_tracker.get_fraction_cov_atleast(2) == approx(11/17)
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr1", "chr2"]) == approx(9/13)
    
    assert cov_tracker.avg_chrom_coverage("chr1") == covs_per_bp_per_chrom["chr1"].mean()
    assert cov_tracker.avg_chrom_coverage("chr1", metric="median") == np.median(covs_per_bp_per_chrom["chr1"])
    
    assert cov_tracker.get_chrom_lens() == {chrom: len(covs) for chrom, covs in covs_per_bp_per_chrom.items()}

def test_covtracker_plotting():
    covs_per_bp_per_chrom = {"chr1": np.array([1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 3, 3, 0, 0, 5, 5, 5, 5, 1, 1, 1, 0, 0, 0, 0, 0]), "chr2": np.array([1, 1, 2, 2, 3, 3, 4, 4, 4])}
    cov_tracker = NanoSimCoverageTracker(covs_per_bp_per_chrom, {chrom: len(covs) for chrom, covs in covs_per_bp_per_chrom.items()}, blocksize=1)
    
    cov_tracker.plot_state("line", target_coverage=4, figsize=(6, 6)).set_layout_engine(TightLayoutEngine())
    cov_tracker.plot_state("hist", target_coverage=4, figsize=(6, 4)).set_layout_engine(TightLayoutEngine())


##### test adding reads #####

def test_covtracker_add_reads():
    covs_per_bp_per_chrom = {
        "chr1": np.array([1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 3, 3, 0, 0, 5, 5, 5, 5, 1, 1, 1, 0, 0, 0, 0, 0]),
        "chr2": np.array([1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 3, 3, 0, 0, 5, 5, 5, 5, 1, 1, 1, 0, 0, 0, 0, 0])
    }
    cov_tracker1 = NanoSimCoverageTracker(copy.deepcopy(covs_per_bp_per_chrom), {chrom: len(covs) for chrom, covs in covs_per_bp_per_chrom.items()}, blocksize=1)
    cov_tracker2 = NanoSimCoverageTracker(copy.deepcopy(covs_per_bp_per_chrom), {chrom: len(covs) for chrom, covs in covs_per_bp_per_chrom.items()}, blocksize=1)
    
    cov_tracker1.add_reads(["chr1_5_aligned_proc0:0_F_0_7_0", "chr2_15_aligned_proc0:0_R_0_7_0"], verbose=False)
    covs_per_bp_per_chrom["chr1"][5:5+7] += 1
    covs_per_bp_per_chrom["chr2"][15:15+7] += 1 # ref position with respect to forward strand
    assert_dict_with_np_arrays_equal(cov_tracker1.coverage_per_chrom, covs_per_bp_per_chrom)
    
    cov_tracker2.add_read("chr1_5_aligned_proc0:0_F_0_7_0")
    cov_tracker2.add_read("chr2_15_aligned_proc0:0_R_0_7_0")
    assert_dict_with_np_arrays_equal(cov_tracker2.coverage_per_chrom, covs_per_bp_per_chrom)

def test_PafCoverageTracker(tmp_path):
    dummy_paf = tmp_path / "dummy.paf"
    with open(dummy_paf, "w") as f:
        f.write(dedent(
        """\
        read1	474	9	466	-	chr1	1000000	275617	276074	432	457	60	tp:A:P	cm:i:77	s1:i:432	s2:i:51	dv:f:0.0009	rl:i:70
        read2	373	3	364	-	chr2	2000000	564935	565296	361	361	60	tp:A:P	cm:i:74	s1:i:361	s2:i:0	dv:f:0.0009	rl:i:0
        read4	795	0	792	+	chr2	2000000	226989	227781	792	792	60	tp:A:P	cm:i:144	s1:i:792	s2:i:0	dv:f:0.0005	rl:i:0"""
        ))
        
    cov_tracker = PafCoverageTracker.empty_from_paf(paf_file=dummy_paf)
    assert cov_tracker.chrom_lens == {"chr1": 1000000, "chr2": 2000000}
    cov_tracker.get_chrom_start_len("read1") == ("chr1", 275617, 457)
    cov_tracker.get_chrom_start_len("read2") == ("chr2", 564935, 361)
    cov_tracker.get_chrom_start_len("read3") is None
    
##### CovTracker with blocks

def test_covblocktracker_from_empty(tmp_path):
    ref_genome_file = tmp_path / "ref_genome.fa"
    with open(ref_genome_file, "w") as f:
        f.write(dedent("""\
            > chr1
            AAATTTTGGGTTTCCCCA
            > chr2
            CCTTAAATTAAAATTCCCCCCCGGGGG
            """))
        
    assert_dict_with_np_arrays_equal(NanoSimCoverageTracker.empty_from_ref_genomes(ref_genome_file, blocksize=6).coverage_per_chrom, {"chr1": np.zeros(3), "chr2": np.zeros(5)}) # 18/6, 27/6
    
    assert_dict_with_np_arrays_equal(NanoSimCoverageTracker.empty_from_lens({"chr1": 18, "chr2": 27}, blocksize=6).coverage_per_chrom, {"chr1": np.zeros(3), "chr2": np.zeros(5)})

def test_covblocktracker_math():
    # test functions that compute something
    covs_per_bp_per_chrom = {
        "chr1": np.array([4, 8, 12, 11, 3, 15, 7, 0, 2]),
        "chr2": np.array([5, 3, 2, 2, 7]),
        "chr3": np.array([1, 1, 3, 3])
    }
    chrom_lens = {"chr1": 25, "chr2": 14, "chr3": 12}
    blocksize = 3
    cov_tracker = NanoSimCoverageTracker(copy.deepcopy(covs_per_bp_per_chrom), chrom_lens, blocksize=blocksize)
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr1"]) == approx((5*blocksize + 1)/25) # last block has size 1
    assert cov_tracker.get_fraction_cov_atleast(3, ["chr1"]) == approx(3*blocksize/25)
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr2"]) == approx(2/14) # last block has size 2
    assert cov_tracker.get_fraction_cov_atleast(3, ["chr3"]) == approx(0/12)
    assert cov_tracker.get_fraction_cov_atleast(2) == approx((5*blocksize + 1 + 2 + 0)/(25+14+12))
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr1", "chr2"]) == approx((5*blocksize + 1 + 2)/(25+14))
    
    assert cov_tracker.avg_chrom_coverage("chr1") == np.array([4/3, 8/3, 12/3, 11/3, 3/3, 15/3, 7/3, 0/3] * 3 + [2]).mean()
    assert cov_tracker.avg_chrom_coverage("chr1", metric="median") == np.median(np.array([4/3, 8/3, 12/3, 11/3, 3/3, 15/3, 7/3, 0/3] * 3 + [2]))
    
    assert cov_tracker.get_chrom_lens() == chrom_lens
    
def test_covblocktracker_math_large():
    # only test it for the covtracker with blocks since the one without blocks uses uint8
    
    # check no overflow occurs with large numbers (np.uint32)
    blocksize = 2**29
    covs_per_bp_per_chrom = {
        "chr1": np.array([5, 2 * blocksize + 1, 3 * blocksize - 1, blocksize, blocksize + 2], dtype=np.uint32),
        "chr2": np.array([3 * blocksize, 2 * blocksize], dtype=np.uint32),
    }
    chrom_lens = {"chr1": 4 * blocksize + blocksize//2, "chr2": blocksize + 10}
    cov_tracker = NanoSimCoverageTracker(copy.deepcopy(covs_per_bp_per_chrom), chrom_lens, blocksize=blocksize)
    assert cov_tracker.avg_chrom_coverage("chr1") == approx(covs_per_bp_per_chrom["chr1"].sum() / chrom_lens["chr1"])
    
    block_lengths = np.array([blocksize, blocksize, blocksize, blocksize, chrom_lens["chr1"] - 4 * blocksize])
    assert cov_tracker.get_fraction_cov_atleast(2, ["chr1"]) == approx(np.average(covs_per_bp_per_chrom["chr1"] / block_lengths >= 2, weights=block_lengths/block_lengths.sum()))
    assert cov_tracker.get_fraction_cov_atleast(3, ["chr1"]) == approx(np.average(covs_per_bp_per_chrom["chr1"] / block_lengths >= 3, weights=block_lengths/block_lengths.sum()))
    block_lengths = np.array([blocksize, blocksize, blocksize, blocksize, chrom_lens["chr1"] - 4 * blocksize, blocksize, chrom_lens["chr2"] - blocksize])
    assert cov_tracker.get_fraction_cov_atleast(2) == approx(np.average(np.concatenate(list(covs_per_bp_per_chrom.values())) / block_lengths >= 2, weights=block_lengths/block_lengths.sum()))
    assert cov_tracker.get_fraction_cov_atleast(3) == approx(np.average(np.concatenate(list(covs_per_bp_per_chrom.values())) / block_lengths >= 3, weights=block_lengths/block_lengths.sum()))
    
    assert cov_tracker.get_chrom_lens() == chrom_lens
    
def test_covblocktracker_plotting():
    covs_per_block_per_chrom = {"chr1": np.array([4, 8, 12, 11, 3, 15, 7, 0, 1])}
    blocksize = 3
    cov_tracker = NanoSimCoverageTracker(covs_per_block_per_chrom, {"chr1": 25}, blocksize=blocksize)
    
    cov_tracker.plot_state("line", target_coverage=4, figsize=(6, 6)).set_layout_engine(TightLayoutEngine())
    cov_tracker.plot_state("hist", target_coverage=4, figsize=(6, 4)).set_layout_engine(TightLayoutEngine())
    