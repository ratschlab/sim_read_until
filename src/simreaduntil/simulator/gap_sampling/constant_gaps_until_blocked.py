"""
Implements constant short and long gaps, and blocks the channel after a certain number of bps have been read.
"""


from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_long_gap_threshold
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler, compute_shortgap_longgap_longgapprob, fit_nb_bps_until_blocked, fit_seq_time_beta, fit_seq_time_mle
from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_median_read_delay


import numpy as np


from typing import Tuple


class ConstantGapsUntilBlocked(GapSampler):
    """
    Chooses short and long gaps of constant length, where a long gap is chosen with some probability, until the channel is blocked.
    
    Called constant_gaps in the paper
    
    Args:
        short_gap_length: length of short gaps
        long_gap_length: length of long gaps
        prob_long_gap: probability of a long gap
        time_until_blocked: number of bps until the channel is blocked (one it is at least this number)
        read_delay: delay between read starting and first bp being read
    """
    def __init__(self, short_gap_length, long_gap_length, prob_long_gap, time_until_blocked, read_delay) -> None:
        super().__init__()

        self.short_gap_length = short_gap_length
        self.long_gap_length = long_gap_length
        self.prob_long_gap = prob_long_gap
        self.time_until_blocked = time_until_blocked
        self.read_delay = read_delay

    def save(self, save_path):
        np.savez_compressed(save_path, short_gap_length=self.short_gap_length, long_gap_length=self.long_gap_length, prob_long_gap=self.prob_long_gap, time_until_blocked=self.time_until_blocked, read_delay=self.read_delay)
        
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None):
        """
        Compute parameters from a sequencing summary dataframe
        Call add_previous_gap_duration before calling this function
        It assumes that each channel becomes broken after its final read.
        Returns:
            function to create a gap sampler, so it is flexible with respect to the number of channels
        """
        if read_delay is None:
            read_delay = compute_median_read_delay(seqsum_df)

        # dist = fit_nb_bps_until_blocked(seqsum_df)
        # dist = fit_seq_time_beta(seqsum_df)
        dist = fit_seq_time_mle(seqsum_df)

        gap_durations = seqsum_df["prev_gap_duration"]
        long_gap_threshold = compute_long_gap_threshold(gap_durations)
        short_gap_median, long_gap_median, prob_long_gap = compute_shortgap_longgap_longgapprob(gap_durations, long_gap_threshold)

        channel_nb = 0
        def get_gap_sampler(random_state=None):
            time_until_blocked = dist.rvs(random_state=random_state)
            nonlocal channel_nb
            channel_nb += 1
            return channel_nb, cls(short_gap_median, long_gap_median, prob_long_gap=prob_long_gap, time_until_blocked=time_until_blocked, read_delay=read_delay)

        return get_gap_sampler

    def __repr__(self) -> str:
        # repr(random_state) is not very informative (does not show seed, so we store it separately and display it here)
        return f"""ConstantGapsUntilBlocked(short_gap={self.short_gap_length}, long_gap={self.long_gap_length}, prob_long_gap={self.prob_long_gap}, time_until_blocked={self.time_until_blocked}, read_delay={self.read_delay})"""

    def _check_sim_params(self):
        assert isinstance(self.short_gap_length, (int, float))
        assert self.short_gap_length >= 0

        assert isinstance(self.long_gap_length, (int, float))
        assert self.long_gap_length >= 0

        assert isinstance(self.prob_long_gap, (int, float))
        assert 1 >= self.prob_long_gap >= 0

        assert isinstance(self.read_delay, (int, float))
        assert self.read_delay >= 0

    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        return self.read_delay

    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        # if channel_stats.time_active_without_mux_scan >= self.time_until_blocked:
        # if channel_stats.reads.number_bps_read >= self.time_until_blocked:
        t = channel_stats.time_active_without_mux_scan
        if t >= self.time_until_blocked:
            return GapSampler.GapType.Broken, np.inf
        elif random_state.uniform() < self.prob_long_gap:
            return GapSampler.GapType.Long, self.long_gap_length
        else:
            return GapSampler.GapType.Short, self.short_gap_length