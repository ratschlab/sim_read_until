import functools
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pytest
from pytest import approx
from simreaduntil.seqsum_tools.seqsum_plotting import plot_channels_over_time
from simreaduntil.seqsum_tools.seqsum_preprocessing import add_previous_gap_duration, sort_and_clean_seqsum_df, compute_median_pore_speed
from simreaduntil.shared_utils.plotting import make_tight_layout
from simreaduntil.shared_utils.utils import subset_dict, get_some_value_from_dict
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked

from simreaduntil.simulator.gap_sampling.gap_sampling import RandomGapSampler, compute_prob_long_gap, compute_shortgap_longgap_longgapprob, get_gaps_df_including_last_gap, restrict_gaps_df_to_window
from simreaduntil.simulator.gap_sampling.inactive_active_gaps_replication import SingleChannelInactiveActiveReplicator, get_read_durations_per_channel, plot_inactive_periods, plot_inactive_periods_single_channel
from simreaduntil.simulator.gap_sampling.gap_sampler_per_window_until_blocked import compute_time_and_aggregation_windows, GapSamplerPerWindowUntilBlocked
from simreaduntil.simulator.gap_sampling.rolling_window_gap_sampler import RollingWindowGapSamplerPerChannel
from simreaduntil.simulator.simfasta_to_seqsum import convert_simfasta_dir_to_seqsum, convert_simfasta_to_seqsum
from simreaduntil.simulator.simulator import assign_read_durations_to_channels, run_simulator_from_sampler_per_channel, run_simulator_from_sampler_per_channel_parallel
from simreaduntil.simulator.simulator_params import SimParams

get_random_seqsum_df = pytest.helpers.get_random_seqsum_df

def test_compute_prob_long_gap():
    assert compute_prob_long_gap(2., 2., 0.3) == approx(0.3) # preserved
    assert compute_prob_long_gap(2., 2., 0.5) == approx(0.5) # preserved
    assert compute_prob_long_gap(2., 3., 1e-16) == approx(0.)
    assert compute_prob_long_gap(2.0, 3.0, 1.) == approx(1.)
    assert compute_prob_long_gap(2.0, 3.0, 0.5) < 0.45 # prob of choosing a long gap should decrease, i.e. smaller 0.5 (test 0.45 here)

    compute_shortgap_longgap_longgapprob(np.array([1, 5, 8, 2, 3, 9]), long_gap_threshold=4)

def test_get_gaps_df_including_last_gap():
    df_expected = pd.DataFrame.from_records(
        [(1, 0., 1.5), (1, 5., 7.), (2, 0., 3.), (2, 7., 9.)],
        columns=["channel", "gap_start", "gap_end"],
    )
    df_expected["gap_duration"] = df_expected["gap_end"] - df_expected["gap_start"]
    df_expected
    pd.testing.assert_frame_equal(
        df_expected,
        get_gaps_df_including_last_gap(
            pd.DataFrame.from_records(
                [(1, 1.5, 5., 1.5), (1, 7., 9., 2.), (2, 3., 7., 3.)],
                columns=["channel", "start_time", "end_time", "prev_gap_duration"],
            )
        ),
        check_like=True
    )

def test_restrict_gaps_df_to_window():
    pd.testing.assert_frame_equal(
        restrict_gaps_df_to_window(
            pd.DataFrame.from_records([(0, 3), (2, 5), (2, 7), (6, 12), (10, 12), (3, 13), (15, 20)], columns=["gap_start", "gap_end"]),
            5, 10
        ).reset_index(drop=True), 
        pd.DataFrame.from_records([(2, 5), (2, 7), (6, 12), (10, 12), (3, 13)], columns=["gap_start", "gap_end"]),
        check_like=True
    )

def run_simulator_with_gap_sampler(gap_sampler_from_seqsum_fcn, df_read, seq_end_time, run_dir, n_channels):
    """
    Run simulator for a set of channels
    
    Args:
        gap_sampler_from_seqsum_fcn (function): function that takes a seqsum_df and returns a function that creates gap samplers
        df_read (pd.DataFrame): seqsum_df as read from the file
        seq_end_time (float): end time of the sequence
        run_dir (pathlib.Path): directory to save the results
        n_channels (int): number of channels to simulate
    """
    seqsum_df = sort_and_clean_seqsum_df(df_read, min_gap_duration=1e-8)
    seqsum_df = add_previous_gap_duration(seqsum_df, seq_start_time=0)

    random_state = np.random.default_rng(3)
    
    gap_sampler_maker = gap_sampler_from_seqsum_fcn(seqsum_df)
    gap_samplers = dict(gap_sampler_maker(random_state=random_state) for _ in range(n_channels))
    read_durations_per_channel = get_read_durations_per_channel(seqsum_df)
    # possibly restrict if simulating less channels
    read_durations_per_channel = assign_read_durations_to_channels(read_durations_per_channel.values(), channel_names=gap_samplers.keys())
    
    if isinstance(gap_samplers[1], SingleChannelInactiveActiveReplicator):
        plot_inactive_periods_single_channel(gap_samplers[2].inactive_active_periods_tracker)
        gap_samplers[1].inactive_active_periods_tracker.get_all_gaps()
        plot_inactive_periods(gap_samplers)

    run_dir.mkdir()
    gap_sampler_dir = run_dir / "gap_samplers"
    gap_sampler_dir.mkdir()
    for (channel, gap_sampler) in gap_samplers.items():
        gap_sampler.save(gap_sampler_dir / f"gap_sampler_{channel}.npz")
        break # just test one gap sampler can be saved
    # check loading works
    type(get_some_value_from_dict(gap_samplers)).load(gap_sampler_dir / f"gap_sampler_{channel}.npz")
    
    mk_run_dir = run_dir / "reads"
    sim_params = SimParams(gap_samplers=gap_samplers, bp_per_second=compute_median_pore_speed(seqsum_df))
    # run_simulator_from_sampler_per_channel(
    #     mk_run_dir, sim_params=sim_params, read_durations_per_channel=read_durations_per_channel, cycle_read_durations=True, seq_end_time=seq_end_time
    # )
    run_simulator_from_sampler_per_channel_parallel(
        mk_run_dir, sim_params=sim_params, read_durations_per_channel=read_durations_per_channel, cycle_read_durations=True, seq_end_time=seq_end_time
    )
    
    sim_seq_summary_file = run_dir / "sequencing_summary.txt"
    convert_simfasta_dir_to_seqsum(mk_run_dir, seqsummary_filename=sim_seq_summary_file)
    sim_seqsum_df = sort_and_clean_seqsum_df(pd.read_csv(sim_seq_summary_file, sep="\t"))
    assert len(sim_seqsum_df) > 0
    sim_seqsum_df = add_previous_gap_duration(sim_seqsum_df, seq_start_time=0)
    ax = plot_channels_over_time(sim_seqsum_df, figsize=(20, 10))
    # save_fig_and_pickle(ax.figure, "simulated_run.png") # for the README.md, markdown does not natively render pdf, keep commented out
    plt.close(ax.figure)

def _get_gap_sampler_per_window_until_blocked(seqsum_df, read_delay=None):
    # to set the random state and time windows
    time_and_aggregation_windows = compute_time_and_aggregation_windows(seq_end_time=seqsum_df["end_time"].max(), nb_windows=10, fraction_overlap=0.5)
    return GapSamplerPerWindowUntilBlocked.from_seqsum_df(seqsum_df, read_delay=read_delay, time_and_aggregation_windows=time_and_aggregation_windows)
    
@pytest.mark.parametrize("gap_sampler_from_seqsum_fcn", [
    functools.partial(RollingWindowGapSamplerPerChannel.from_seqsum_df, n_channels_full=10), # contains 9 and 6 channels respectively
    RandomGapSampler.from_seqsum_df, 
    ConstantGapsUntilBlocked.from_seqsum_df,
    _get_gap_sampler_per_window_until_blocked,
    SingleChannelInactiveActiveReplicator.from_seqsum_df
])
def test_run_replication(shared_datadir, tmp_path, gap_sampler_from_seqsum_fcn):
    # test that parameters can be learnt from a run and then be simulated
    
    # todo2: np raises a DeprecationWarning, but cannot see traceback
    # np.seterr(all='raise') # todo
    # import warnings
    # warnings.filterwarnings("error")
    
    # note: this results in a fiterror, probably because there is too little data; it works on full sequencing summary files though
    sequencing_summary_file = shared_datadir / "zymo_short_seqsum.txt"
    df_read = pd.read_csv(sequencing_summary_file, sep="\t")
    run_simulator_with_gap_sampler(gap_sampler_from_seqsum_fcn, df_read=df_read, seq_end_time=15_000, run_dir=tmp_path / "run_replication_real", n_channels=6)
    
    df_read = get_random_seqsum_df()
    run_simulator_with_gap_sampler(gap_sampler_from_seqsum_fcn, df_read=df_read, seq_end_time=15_000, run_dir=tmp_path / "run_replication_sim", n_channels=6)
    