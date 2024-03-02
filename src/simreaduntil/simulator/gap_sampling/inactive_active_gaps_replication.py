"""
Replicate a real run by extracting active and inactive periods, as well as read times.
"""

import itertools
import pickle
from typing import Dict, Generator, List, Tuple, Any
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

import numpy as np
import tqdm
from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_long_gap_threshold, get_gaps_single_channel
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.plotting import make_tight_layout
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.gap_sampling import CompletelyBrokenGapSampler, GapSampler
from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_median_read_delay

logger = setup_logger_simple(__name__)
"""module logger"""

def _interleave(a, b):
    """
    Interleave two arrays where len(a) = len(b) or len(a) = len(b) + 1
    """
    assert len(a) in [len(b), len(b) + 1]
    if len(b) == 0:
        res = np.array([])
    else:
        res = np.concatenate(list(zip(a, b)))
    if len(a) == len(b):
        return res
    return np.concatenate((res, [a[-1]]))

class ChannelInactiveActivePeriodsTracker:
    """
    Keep track of active and inactive periods for a single channel, similar to UNCALLED
    
    The channel consists of alternating active and inactive periods. Each inactive period is also called a long gap.
    Each active period has a duration and consists of short gaps that are added between the reads.
    """
    def __init__(self):
        self.inactive_period_lengths = []
        self.short_gaps_per_active_period = []
        self.active_period_lengths = []
        self.first_period_is_active = None # will be set on first call to add_inactive_period or add_active_period
    
    def add_inactive_period(self, long_gap):
        """
        Add an inactive period
        """
        if self.first_period_is_active is None:
            self.first_period_is_active = False
            
        assert len(self.inactive_period_lengths) == len(self.short_gaps_per_active_period) - int(self.first_period_is_active), "must add short and long gaps in alternation"
        self.inactive_period_lengths.append(long_gap)
        
    def add_active_period(self, active_period_length, short_gaps):
        """
        Add an active period
        
        An active period may also have no short gaps!
        """
        if self.first_period_is_active is None:
            self.first_period_is_active = True
            
        # assert len(short_gaps) > 0, "active period must have at least one short gap"
            
        assert len(self.short_gaps_per_active_period) == len(self.inactive_period_lengths) - int(not self.first_period_is_active), "must add short and long gaps in alternation"
        self.short_gaps_per_active_period.append(short_gaps)
        self.active_period_lengths.append(active_period_length)
    
    @property    
    def num_periods(self):
        return len(self.inactive_period_lengths) + len(self.active_period_lengths)
    
    @property
    def full_duration(self):
        """
        Full duration of that can be simulated
        
        This is a lower bound because a period may run longer due to an active read.
        """
        return sum(self.inactive_period_lengths) + sum(self.active_period_lengths)
    
    def check_valid_period(self, period_idx):
        """Check that period_idx is in range"""
        assert 0 <= period_idx < self.num_periods, f"period_idx={period_idx} is out-of-range [0, {len(self.inactive_period_lengths) + len(self.active_period_lengths)})"
        
    def is_valid_period(self, period_idx):
        """Whether the period_idx is in range"""
        return 0 <= period_idx < self.num_periods
    
    def is_active_period(self, period_idx):
        """
        Whether period is active
        
        The channel is a succession of either:
        - active, inactive, active, inactive, ...
        - inactive, active, inactive, active, ...
        """
        return period_idx % 2 != self.first_period_is_active
    
    def is_last_period(self, period_idx):
        """Whether this is the last period"""
        return period_idx == self.num_periods - 1
    
    def get_period_duration(self, period_idx):
        """
        Get length of ith period
        """
        if self.is_active_period(period_idx):
            return self.active_period_lengths[period_idx // 2]
        else:
            return self.inactive_period_lengths[period_idx // 2]
        
    def get_inactive_period_positions(self):
        """
        Get inactive periods positions [(start, end), ]
        """
        if self.first_period_is_active:
            durations = _interleave(self.active_period_lengths, self.inactive_period_lengths)
        else:
            durations = _interleave(self.inactive_period_lengths, self.active_period_lengths)
        period_positions = np.concatenate(([0], durations.cumsum())) # (start1, end1, start2, end2, ...)
        
        period_positions = period_positions[int(self.first_period_is_active):] # get rid of first period if it is active
        return period_positions[:2 * (len(period_positions)//2)].reshape(-1, 2) # get rid of last period if it is active (i.e. make even length), then reshape
    
    def get_all_gaps(self):
        # returns all gaps in the order they appear
        return np.concatenate([
            np.concatenate((x, y)) if self.first_period_is_active else np.concatenate((y, x))
            for (x, y) in itertools.zip_longest(self.short_gaps_per_active_period, [np.array([x]) for x in self.inactive_period_lengths], fillvalue=np.array([]))
        ])
        
    def __repr__(self):
        return f"GapsBetweenLong({len(self.inactive_period_lengths)} long_gaps, {len(self.short_gaps_per_active_period)} short_gaps_per_active, first_gap_is_active={self.first_period_is_active})"
    
def extract_active_inactive_periods(df_single, seq_start_time, seq_end_time, long_gap_threshold=None) -> ChannelInactiveActivePeriodsTracker:
    """
    Extract active and inactive periods for a single channel
    
    Args:
        df_single: dataframe with reads for a single channel (start_time and end_time columns)
        seq_start_time: start time of the whole sequencing run
        seq_end_time: end time of the whole sequencing run
        long_gap_threshold: threshold for long gaps (if None, use median + std)
        
    Returns:
        object containing the lengths of active and inactive periods, as well as short gaps within each period
    """
    gap_durations, gap_starts = get_gaps_single_channel(df_single, seq_start_time=seq_start_time)
    gap_ends = gap_starts + gap_durations

    if long_gap_threshold is None:
        long_gap_threshold = compute_long_gap_threshold(gap_durations)
    
    channel_gaps_tracker = ChannelInactiveActivePeriodsTracker()
    long_gap_indices = np.where(gap_durations > long_gap_threshold)[0]
        
    # long gap is an inactive period, in between are active periods with short gaps
    # channel always starts with a gap (either short or long). If it is long, it starts with an inactive period, otherwise with an active period.
    
    period_start = seq_start_time # start of last active/inactive period
    short_gap_idx_start = 0
        
    for next_long_gap_idx in long_gap_indices:
        if next_long_gap_idx > 0:
            # first gap is short, so we start with an active period
            if short_gap_idx_start < next_long_gap_idx:
                short_gaps = gap_durations[short_gap_idx_start:next_long_gap_idx]
            else:
                # logger.warning(f"Two long gaps in a row with no short gaps in between at index {next_long_gap_idx}, using a default value")
                short_gaps = []
            channel_gaps_tracker.add_active_period(gap_starts[next_long_gap_idx] - period_start, short_gaps)
        
        channel_gaps_tracker.add_inactive_period(gap_durations[next_long_gap_idx])
        
        period_start = gap_ends[next_long_gap_idx]
        short_gap_idx_start = next_long_gap_idx + 1
    
    # add short gaps in last active period, if any
    short_gaps = gap_durations[short_gap_idx_start:]
    if len(short_gaps) > 0:
        channel_gaps_tracker.add_active_period(seq_end_time - period_start, short_gaps)
    
    return channel_gaps_tracker

def plot_inactive_periods_single_channel(inactive_active_periods_tracker: ChannelInactiveActivePeriodsTracker, y_pos=0):
    """
    Plot inactive periods for a single channel
    """
    fig, ax = plt.subplots()
    
    segments =[((start, y_pos), (end, y_pos)) for (start, end) in inactive_active_periods_tracker.get_inactive_period_positions()]
    line_collection = LineCollection(segments, linewidths=2, colors="blue")
    ax.add_collection(line_collection)
    ax.autoscale()


def get_read_durations_per_channel(seqsum_df) -> Dict[Any, np.ndarray]:
    """Return a dictionary with read durations per channel"""
    return {channel: df_single["duration"].values for (channel, df_single) in seqsum_df.groupby("channel", observed=True)}


class SingleChannelInactiveActiveReplicator(GapSampler):
    """
    Replicate a single channel by replicating the lengths of active and inactive periods as much as possible.
    
    Whenever the length of an active period is exceeded, move to the next period.
    Whenever an inactive period finishes, you must call mark_long_gap_end. This moves to the next (active) period.
    Within an active period, the short gaps are recycled from the observed short gaps within the active period.
    It ignores the time spent in mux scans.
    
    Called gap_replication in the paper
    """
    def __init__(self, inactive_active_periods_tracker: ChannelInactiveActivePeriodsTracker, read_delay) -> None:
        super().__init__()

        self.inactive_active_periods_tracker = inactive_active_periods_tracker
        self.read_delay = read_delay
        if len(self.inactive_active_periods_tracker.short_gaps_per_active_period) == 0:
            # no active periods
            self.default_short_gap = 0
        else:
            # gap when an active period has no short gaps, no short gaps -> NaN, replace by 0
            self.default_short_gap = np.nan_to_num(np.median(np.concatenate(self.inactive_active_periods_tracker.short_gaps_per_active_period)), 0)
        assert self.default_short_gap >= 0, f"Encountered negative default gap {self.default_short_gap}"

        self.cur_period_idx = -1
        self.cur_period_end = 0
        self._move_to_next_period()

    def save(self, save_path):
        # np.savez_compressed(save_path, inactive_active_periods_tracker=self.inactive_active_periods_tracker.dict_for_np_saving(), read_delay=self.read_delay)
        # use pickle since the nested object is not a numpy array, so would need methods to load and save itself
        with open(save_path, "wb") as f:
            pickle.dump({"read_delay": self.read_delay, "inactive_active_periods_tracker": self.inactive_active_periods_tracker}, f)
        
    @classmethod
    def load(cls, save_path):
        with open(save_path, "rb") as f:
            return cls(**pickle.load(f))
        
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None):
        """
        Wrapper which returns the gap sampler
        
        Calling this function n_channel times returns gap samplers for all channels
        Calling it more often produces an error.
        
        Args:
            seqsum_df: (preprocessed) dataframe from a sequencing summary file
            read_delay: read delay of each read; if None, it is computed from the sequencing summary
        
        Returns:
            function to create a gap sampler
        """

        def _extract_gap_samplers(seqsum_df, read_delay) -> Generator[Tuple[Any, GapSampler], None, None]:
            """
            Extract gap patterns from a sequencing summary to replicate the corresponding run
            See SingleChannelInactiveActiveReplicator
            Note: If a channel has no reads, it is skipped.
            
            Once the number of active channels is depleted, it returns a CompletelyBrokenGapSampler.
            
            Args:
                seqsum_df: (preprocessed) dataframe from a sequencing summary file
                read_delay: read delay of each read
            """
            seq_start_time = 0
            seq_end_time = seqsum_df["end_time"].max()
            
            df_channel_groups = seqsum_df.groupby("channel", observed=True)

            # global long gap threshold because this gives more data (a single channel may have very few gaps!)
            long_gap_threshold = compute_long_gap_threshold(seqsum_df["prev_gap_duration"])

            channel_nb = 0
            for (channel, df_single) in tqdm.tqdm(df_channel_groups, desc="Extracting replication parameters from channels"):
                channel_gaps_tracker = extract_active_inactive_periods(df_single, seq_start_time=seq_start_time, seq_end_time=seq_end_time, long_gap_threshold=long_gap_threshold)

                channel_nb += 1
                yield (channel_nb, SingleChannelInactiveActiveReplicator(channel_gaps_tracker, read_delay=read_delay))
                
            while True:
                channel_nb += 1
                yield (channel_nb, CompletelyBrokenGapSampler())

        if read_delay is None:
            read_delay = compute_median_read_delay(seqsum_df)
        gap_samplers_iter = iter(_extract_gap_samplers(seqsum_df, read_delay=read_delay))
        def wrapper(random_state=None):
            return next(gap_samplers_iter)
        return wrapper

    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        return self.read_delay

    def mark_long_gap_end(self, channel_stats: ChannelStats):
        # channel_stats: with stats at the time the long gap ends
        assert not self.inactive_active_periods_tracker.is_active_period(self.cur_period_idx)
        self.cur_period_end = channel_stats.time_active_without_mux_scan
        self._move_to_next_period()

    def _move_to_next_period(self):
        """
        Move to the next period
        """
        self.cur_period_idx += 1
        if self.inactive_active_periods_tracker.is_valid_period(self.cur_period_idx):
            self.cur_period_end += self.inactive_active_periods_tracker.get_period_duration(self.cur_period_idx)
        else:
            # channel broken
            self.cur_period_end = np.inf
        self.active_period_short_gap_idx = -1

    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        """
        Sample the next gap
        
        Once a long gap ends, you should call mark_long_gap_end.
        This function is usually called after a read
        gap - read - gap - read ...
        If the current period is over, it moves to the next period.
        When the last period is over, it is followed by a ChannelBroken gap.
        """

        read_end_time = channel_stats.time_active_without_mux_scan

        if self.inactive_active_periods_tracker.is_active_period(self.cur_period_idx):
            # if active period is over, move to next period
            if read_end_time >= self.cur_period_end:
                self.cur_period_end = read_end_time
                self._move_to_next_period()
        else:
            assert read_end_time < self.cur_period_end, "inactive period should be ended with mark_long_gap_end"

        if not self.inactive_active_periods_tracker.is_valid_period(self.cur_period_idx):
            # no more periods left, just go into ChannelBroken state
            return GapSampler.GapType.Broken, np.inf
        elif self.inactive_active_periods_tracker.is_active_period(self.cur_period_idx):
            self.active_period_short_gap_idx += 1
            short_gaps = self.inactive_active_periods_tracker.short_gaps_per_active_period[self.cur_period_idx // 2]
            short_gap = short_gaps[self.active_period_short_gap_idx % len(short_gaps)] if len(short_gaps) > 0 else self.default_short_gap
            return GapSampler.GapType.Short, short_gap
        else:
            return GapSampler.GapType.Long, self.inactive_active_periods_tracker.get_period_duration(self.cur_period_idx)

    def is_broken(self):
        # whether the channel is in a broken state (no more periods left)
        return not self.inactive_active_periods_tracker.is_valid_period(self.cur_period_idx)

    @property
    def full_duration(self):
        return self.inactive_active_periods_tracker.full_duration

def plot_inactive_periods(gap_samplers: List[SingleChannelInactiveActiveReplicator]):
    """
    Plot active/inactive periods
    """
    fig, ax = plt.subplots(figsize=(20, 10))

    for (channel, gap_sampler) in gap_samplers.items():
        y_pos = channel
        segments =[((start, y_pos), (end, y_pos)) for (start, end) in gap_sampler.inactive_active_periods_tracker.get_inactive_period_positions()]
        line_collection = LineCollection(segments, linewidths=2, colors="blue")
        ax.add_collection(line_collection)

    ax.autoscale()

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Channel")
    ax.set_title(f"Inactive periods over time ({len(gap_samplers)} active channels)")
    make_tight_layout(fig)
