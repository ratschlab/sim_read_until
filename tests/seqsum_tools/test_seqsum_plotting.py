
from textwrap import dedent
from unittest import mock
import warnings
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pytest
from pytest import approx

from simreaduntil.seqsum_tools.seqsum_plotting import add_group_and_reflen_from_nanosim_id, add_group_and_reflen_from_paf, compute_coverage_per_group_df, create_plots_for_seqsum, keep_largest_gaps_only, plot_channels_over_time, plot_coverage_per_group, plot_fraction_states_per_channel, plot_number_channels_per_group_over_time, plot_processed_seqsum, seqsum_add_cols_for_plotting_selseq_performance, set_plots_groupby_column
from simreaduntil.seqsum_tools.seqsum_preprocessing import add_previous_gap_duration
from simreaduntil.seqsum_tools.coverage_tracker import NanoSimCoverageTracker, PafCoverageTracker
from simreaduntil.shared_utils.plotting import ignore_tight_layout_warning

set_plots_groupby_column("unit")

ignore_tight_layout_warning()

get_random_seqsum_df = pytest.helpers.get_random_seqsum_df

def test_add_group_and_reflen_from_nanosim_id():
    seqsum_df = pd.DataFrame({"read_id": ["chr1_10_perfect_proc0:1_F_0_840_0", "chr2_120_perfect_proc0:1_F_0_130_0"]})
    expected_seqsum_df = seqsum_df.copy()
    expected_seqsum_df["chrom"] = ["chr1", "chr2"]
    expected_seqsum_df["nb_ref_bps"] = [840, 130]
    add_group_and_reflen_from_nanosim_id(seqsum_df, group_column="chrom")
    pd.testing.assert_frame_equal(seqsum_df, expected_seqsum_df, check_like=True)

def test_add_group_and_reflen_from_paf(tmp_path):
    seqsum_df = pd.DataFrame({"read_id": ["read1", "read2", "read3"], "sequence_length_template": [400, 300, 210]})
    
    dummy_paf = tmp_path / "dummy.paf"
    with open(dummy_paf, "w") as f:
        f.write(dedent(
        """\
        read1	474	9	466	-	chr1	1000000	275617	276074	432	457	60	tp:A:P	cm:i:77	s1:i:432	s2:i:51	dv:f:0.0009	rl:i:70
        read2	373	3	364	-	chr2	1000000	564935	565296	361	361	60	tp:A:P	cm:i:74	s1:i:361	s2:i:0	dv:f:0.0009	rl:i:0
        read4	795	0	792	+	chr2	1000000	226989	227781	792	792	60	tp:A:P	cm:i:144	s1:i:792	s2:i:0	dv:f:0.0005	rl:i:0"""
        ))
    seqsum_df = add_group_and_reflen_from_paf(seqsum_df, dummy_paf, group_column="chrom")
    expected_seqsum_df = seqsum_df.copy()
    expected_seqsum_df["chrom"] = ["chr1", "chr2", "unmapped"]
    expected_seqsum_df["nb_ref_bps"] = [457, 361, 210]
    pd.testing.assert_frame_equal(seqsum_df, expected_seqsum_df, check_like=True)
    
def test_keep_largest_gaps_only():
    segments_x = [(1, 2), (12, 15), (15.5, 30), (35, 50)]
    # gaps are 10, 0.5, 5, so not sorted
    y_pos = -3
    segments = np.array([[(x_start, y_pos), (x_end, y_pos)] for (x_start, x_end) in segments_x])

    assert keep_largest_gaps_only(segments, 0) == approx(np.array([[segments[0][0], segments[-1][1]]]))
    assert keep_largest_gaps_only(segments, 1) == approx(np.array([[[1, -3], [2, -3]], [[12, -3], [50, -3]]]))
    assert keep_largest_gaps_only(segments, 2) == approx(np.array([[[1, -3], [2, -3]], [[12, -3], [30, -3]], [[35, -3], [50, -3]]]))
    assert keep_largest_gaps_only(segments, 5) == approx(segments)
    
def test_plot_channels_over_time(shared_datadir):
    # test that it is reasonably fast
    seqsummary_filename = shared_datadir / "sim_sequencing_summary.txt"
    seqsum_df = pd.read_csv(seqsummary_filename, sep="\t")
    seqsum_df["end_time"] = seqsum_df["start_time"] + seqsum_df["duration"]
    seqsum_df = add_previous_gap_duration(seqsum_df, seq_start_time=0)
    ax = plot_channels_over_time(seqsum_df, max_num_gaps_per_channel=40)
    plt.close(ax.figure)
    
def test_plot_number_channels_per_group_over_time(shared_datadir):
    seqsummary_filename = shared_datadir / "sim_sequencing_summary.txt"
    seqsum_df = pd.read_csv(seqsummary_filename, sep="\t")
    seqsum_df["end_time"] = seqsum_df["start_time"] + seqsum_df["duration"]
    seqsum_df["unit"] = np.random.default_rng(2).choice(["chr1", "chr2"], p=[0.8, 0.2], size=len(seqsum_df))
    fig = plot_number_channels_per_group_over_time(seqsum_df)
    plt.close(fig)
    
def test_plot_fraction_states_per_channel():
    seqsum_df = pd.DataFrame.from_records([
        (1, 2, 3), (1, 5, 10), (1, 12, 20),
        (2, 12, 15), (2, 19, 23), (2, 30, 40)
    ], columns=["channel", "start_time", "end_time"])
    seqsum_df["duration"] = seqsum_df["end_time"] - seqsum_df["start_time"]
    seqsum_df = add_previous_gap_duration(seqsum_df)

    # channel2: 17, 4, 19; channel1: 14, 6, 20
    ax = plot_fraction_states_per_channel(seqsum_df, long_gap_threshold=5)
    plt.close(ax.figure)
    
@pytest.mark.parametrize("seqsummary_filename", ["sim_sequencing_summary.txt", "zymo_short_seqsum.txt"])
def test_seqsum_plots_nanosim(shared_datadir, tmp_path, seqsummary_filename):
    seqsummary_filename = shared_datadir / seqsummary_filename
    
    cov_thresholds = [1, 2]
    
    save_dir = tmp_path / "seqsum_plots"
    save_dir.mkdir()
    seqsum_df = pd.read_csv(seqsummary_filename, sep="\t")
    
    seqsum_df["nb_ref_bps"] = seqsum_df["sequence_length_template"]
    if "unit" not in seqsum_df:
        seqsum_df["unit"] = np.random.default_rng(2).choice(["chr1", "chr2"], size=len(seqsum_df))
    
    seqsum_df = seqsum_add_cols_for_plotting_selseq_performance(seqsum_df)
    # check that all columns starting with "cum_" are increasing, per group
    assert (seqsum_df.groupby("unit", observed=True)[[col for col in seqsum_df.columns if col.startswith("cum_")]].diff().dropna() >= 0).all(axis=None)
    
    seqsum_df = add_previous_gap_duration(seqsum_df, seq_start_time=0)
    warnings.filterwarnings("ignore", message="Only plotting full reads")
    plot_processed_seqsum(seqsum_df, save_dir=save_dir)
    assert len(list(save_dir.iterdir())) > 0

    if "zymo_short_seqsum.txt" != seqsummary_filename.name:
        # only for NanoSim ids
        for cov_every in [1, 10]:
            save_dir = tmp_path / f"seqsum_plots_cov_every_{cov_every}"
            save_dir.mkdir()
            
            cov_tracker = NanoSimCoverageTracker.empty_from_lens({"chr1": 1_000_000, "chr2": 1_000_000})
            cov_df = compute_coverage_per_group_df(seqsum_df[:1000], cov_tracker=cov_tracker, cov_thresholds=cov_thresholds, coverage_every=cov_every, chrom_column="unit")
            plot_coverage_per_group(cov_df, cov_thresholds=cov_thresholds, save_dir=save_dir)
            
            assert len(list(save_dir.iterdir())) > 0
        
def test_create_plots_for_seqsum_nanosimids(tmp_path, mocker):
    mocker.patch("simreaduntil.seqsum_tools.coverage_tracker.NanoSimCoverageTracker.empty_from_ref_genomes", return_value=NanoSimCoverageTracker.empty_from_lens({"chr1": 1_000_000, "chr2": 1_000_000}))
    
    save_dir = tmp_path / "seqsum_plots"
    save_dir.mkdir()
    seqsum_df, cov_df = create_plots_for_seqsum(get_random_seqsum_df(), nrows=100, ref_genome_path="dummy_ref_path", cov_every=10, cov_thresholds=[1, 2], save_dir=save_dir)
    seqsum_df.sort_values("end_time", inplace=True)
    assert (seqsum_df.groupby("chrom", observed=True)[[col for col in seqsum_df.columns if col.startswith("cum_")]].diff().dropna() >= 0).all(axis=None)
    
def test_create_plots_for_seqsum_nanosimids(tmp_path, mocker):
    # using NanoSim ids
    mocker.patch("simreaduntil.seqsum_tools.coverage_tracker.NanoSimCoverageTracker.empty_from_ref_genomes", return_value=NanoSimCoverageTracker.empty_from_lens({"chr1": 1_000_000, "chr2": 1_000_000}))
    
    save_dir = tmp_path / "seqsum_plots"
    save_dir.mkdir()
    seqsum_df, cov_df = create_plots_for_seqsum(get_random_seqsum_df(), nrows=100, ref_genome_path="dummy_ref_path", group_to_units={"chr1-enriched": ["chr1"]}, cov_every=10, cov_thresholds=[1, 2], save_dir=save_dir)
    seqsum_df.sort_values("end_time", inplace=True)
    assert (seqsum_df.groupby("chrom", observed=True)[[col for col in seqsum_df.columns if col.startswith("cum_")]].diff().dropna() >= 0).all(axis=None)

def test_create_plots_for_seqsum_paf(tmp_path):
    # gets chrom, position, len from paf file rather than by parsing read ids
    
    seqsum_df = get_random_seqsum_df()

    # modify seqsum_df so that a paf is added
    seqsum_df["read_id"] = [f"read_{i+1}" for i in range(len(seqsum_df))]
    chroms = seqsum_df["chrom"]
    del seqsum_df["chrom"]
    del seqsum_df["nb_ref_bps"]
    
    dummy_paf_filename = tmp_path / "dummy.paf"
    with open(dummy_paf_filename, "w") as f:
        for (i, chrom) in enumerate(chroms):
            if i <= 5:
                # skip first few entries
                continue
            start_pos = i * 100
            end_pos = start_pos + 40
            f.write(f"""read_{i+1}	474	9	466	-	{chrom}	1000000	{start_pos}	{end_pos}	{end_pos - start_pos}	457	60	tp:A:P	cm:i:77	s1:i:432	s2:i:51	dv:f:0.0009	rl:i:70""")

    # easier to write to a file than mocking because some calls should call original read_csv function
    seqsum_filename = tmp_path / "seqsum_df.txt"
    seqsum_df.to_csv(seqsum_filename, sep="\t", index=False)
    
    save_dir = tmp_path / "seqsum_plots"
    save_dir.mkdir()
    seqsum_df, cov_df = create_plots_for_seqsum(seqsum_filename, nrows=100, paf_file=dummy_paf_filename, cov_every=10, cov_thresholds=[1, 2], save_dir=save_dir)
    