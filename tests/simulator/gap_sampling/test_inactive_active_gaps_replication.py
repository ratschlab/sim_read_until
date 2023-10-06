
import pytest
from pytest import approx
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from simreaduntil.seqsum_tools.seqsum_preprocessing import get_gaps_single_channel
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler
from simreaduntil.simulator.gap_sampling.inactive_active_gaps_replication import ChannelInactiveActivePeriodsTracker, SingleChannelInactiveActiveReplicator, extract_active_inactive_periods, _interleave

def test_interleave():
    assert all(_interleave([], []) == np.array([]))
    assert all(_interleave([0], []) == np.array([0]))
    assert all(_interleave([0, 2], [1, 3]) == np.array([0, 1, 2, 3]))
    assert all(_interleave([0, 2, 4], [1, 3]) == np.array([0, 1, 2, 3, 4]))
    with pytest.raises(AssertionError):
        _interleave([0, 2], [1, 3, 4])

def test_ChannelInactiveActivePeriods():
    inactive_active_periods = ChannelInactiveActivePeriodsTracker()
    
    inactive_active_periods.add_inactive_period(5)
    
    with pytest.raises(AssertionError):
        inactive_active_periods.add_inactive_period(4)
    
    inactive_active_periods.add_active_period(2, [0.5, 0.3])
    
    with pytest.raises(AssertionError):
        inactive_active_periods.add_active_period(4, [0.5, 0.3])
    
    inactive_active_periods.add_inactive_period(4)
    
    assert inactive_active_periods.num_periods == 3
    assert not inactive_active_periods.first_period_is_active
    assert not inactive_active_periods.is_active_period(0)
    assert inactive_active_periods.is_active_period(1)
    
    assert inactive_active_periods.get_period_duration(1) == 2
    assert inactive_active_periods.get_period_duration(2) == 4
    
    assert inactive_active_periods.get_all_gaps() == approx([5, 0.5, 0.3, 4])
    assert (inactive_active_periods.get_inactive_period_positions() == np.array([[0, 5], [7, 11]])).all()
    
    # inactive, active, inactive: above
    
    # inactive, active
    inactive_active_periods = ChannelInactiveActivePeriodsTracker()
    inactive_active_periods.add_inactive_period(5)
    inactive_active_periods.add_active_period(2, [0.2, 0.3])
    assert not inactive_active_periods.first_period_is_active
    assert (inactive_active_periods.get_inactive_period_positions() == np.array([[0, 5]])).all()
    
    # active, inactive, active
    inactive_active_periods = ChannelInactiveActivePeriodsTracker()
    inactive_active_periods.add_active_period(5, [0.2, 0.3])
    inactive_active_periods.add_inactive_period(2)
    inactive_active_periods.add_active_period(4, [0.2, 0.3])
    assert inactive_active_periods.first_period_is_active
    assert (inactive_active_periods.get_inactive_period_positions() == np.array([[5, 7]])).all()
    
    # active, inactive
    inactive_active_periods = ChannelInactiveActivePeriodsTracker()
    inactive_active_periods.add_active_period(5, [0.2, 0.3])
    inactive_active_periods.add_inactive_period(2)
    inactive_active_periods.add_active_period(4, [0.2, 0.3])
    assert inactive_active_periods.first_period_is_active
    assert (inactive_active_periods.get_inactive_period_positions() == np.array([[5, 7]])).all()

def test_simple_extract_active_inactive_periods():
    # simple toy case
    
    # start with a long gap, end with a short gap
    df_single = pd.DataFrame([
        (6, 9), (9.5, 10.5), (12.5, 13.5), (20.5, 23), (30, 33), (36, 39)
    ], columns=["start_time", "end_time"])
    df_single["start_time"] += 20
    df_single["end_time"] += 20
    df_single["channel"] = 1

    seq_start_time = 20
    seq_end_time = 60
    channel_gaps_tracker = extract_active_inactive_periods(df_single, seq_start_time=seq_start_time, seq_end_time=seq_end_time, long_gap_threshold=5)
    assert len(channel_gaps_tracker.get_all_gaps()) == len(df_single)
    assert channel_gaps_tracker.get_all_gaps() == approx(get_gaps_single_channel(df_single, seq_start_time=seq_start_time)[0])
    
# random run
def test_random_extract_active_inactive_periods():
    random_state = np.random.default_rng(2)
    
    # generate (read, gap, read, gap, read, ..., gap, read)
    gap_durations = random_state.choice([1.0, 2.0, 5.0, 6.0], size=100, p=[0.7/2, 0.7/2, 0.3/2, 0.3/2])
    read_durations = random_state.choice([10, 20], size=len(gap_durations), p=[0.5, 0.5])
    read_starts = np.concatenate(([0], (read_durations + gap_durations).cumsum()[:-1]))
    # add gaps around
    df_single = pd.DataFrame({
        "start_time": read_starts,
        "end_time": read_starts + read_durations,
    })
    df_single["start_time"] += 20
    df_single["end_time"] += 20
    df_single["channel"] = 1
    long_gap_threshold = 4

    # check when the channel starts or ends with a short or long gap
    for (offset_low, offset_high) in [(long_gap_threshold/2, long_gap_threshold/2), (long_gap_threshold/2, 2*long_gap_threshold), (2*long_gap_threshold, long_gap_threshold/2), (2*long_gap_threshold, 2*long_gap_threshold)]:
        seq_start_time = df_single["start_time"].min() - offset_low
        seq_end_time = df_single["end_time"].max() + offset_high
        
        channel_gaps_tracker = extract_active_inactive_periods(df_single, seq_start_time=seq_start_time, seq_end_time=seq_end_time, long_gap_threshold=long_gap_threshold)
        assert len(channel_gaps_tracker.get_all_gaps()) == len(df_single)
        
        gap_durations = get_gaps_single_channel(df_single, seq_start_time=seq_start_time)[0]
        assert channel_gaps_tracker.get_all_gaps() == approx(gap_durations)
        assert len(channel_gaps_tracker.inactive_period_lengths) == sum(gap_durations > long_gap_threshold)
        
        long_gap_indices = np.where(gap_durations > long_gap_threshold)[0]
        assert gap_durations[long_gap_indices] == approx(channel_gaps_tracker.inactive_period_lengths)

def test_channel_active_inactive_replication():
    # test that a channel can be replicated from an existing run
    
    df_single = pd.DataFrame([
        (10, 19), (19.5, 20.5), (22.5, 23.5), (30.5, 33), (40, 42), (55, 58.5)
    ], columns=["start_time", "end_time"])
    df_single["channel"] = 1
    
    seq_start_time = 0
    seq_end_time = 60
    channel_gaps_tracker = extract_active_inactive_periods(df_single, seq_start_time=seq_start_time, seq_end_time=seq_end_time, long_gap_threshold=5)
    
    gaps_replicator = SingleChannelInactiveActiveReplicator(channel_gaps_tracker, read_delay=0)
    channel_stats = ChannelStats(n_channels=1)

    # lengths of active periods are not necessarily same as those in the real run due to different read durations
    gap_durations = []
    while not gaps_replicator.is_broken():
        gap_type, gap_duration = gaps_replicator.sample_next_gap(channel_stats, random_state=None) # raises if no more gap available
        gap_durations.append(gap_duration)
        if gap_type == GapSampler.GapType.Short:
            channel_stats.short_gaps.add_full(gap_duration)
        elif gap_type == GapSampler.GapType.Long:
            channel_stats.long_gaps.add_full(gap_duration)
            gaps_replicator.mark_long_gap_end(channel_stats)
            print(f"Long gap: {gap_duration}")
        else:
            print("Channel broken")
        channel_stats.reads.add_full(2.4, nb_new_bps=11, stopped_receiving=False)
        print(f"Finished reads: {channel_stats.reads.finished_number}, time: {channel_stats.time_active}")
            
    default_gap = np.median([0.5, 2.0, 0.5, 2.0, 0.5, 2.0])
    # read duration 2.4, so will recycle gaps
    assert gap_durations == approx([
        10.0, 
        # active period: 2.4 + (0.5 + 2.4) + (2.0 + 2.4) + (0.5 + 2.4) + (2.0 + 2.4) = 17 >= 13.5
        0.5, 2.0, 0.5, 2.0,
        7.0, 
        # active period: 2.4 + (default_gap + 2.4) >= 2.5
        default_gap,
        7.0, 
        # active period: 2.4 >= 2, so no short gap
        13.0,
    ])
    assert channel_stats.reads.finished_number == 5 + 2 + 1 + 1
