"""
Gap sampler that samples from a rolling window, either per channel or with all channels mixed.
"""

from typing import Tuple

import numpy as np
import pandas as pd
from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_long_gap_threshold, compute_median_read_delay
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.simulator.channel_stats import ChannelStats, ElementStats
from simreaduntil.simulator.gap_sampling.gap_sampling import CompletelyBrokenGapSampler, GapSampler, fit_seq_time_mle, get_gaps_df_including_last_gap

logger = setup_logger_simple(__name__)
"""module logger"""

class RollingWindowGapSampler(GapSampler):
    """
    
    Args:
        gaps_df: pd.DataFrame with columns "gap_start", "gap_end", "gap_duration", "gap_type" (long or short)
        window_width: width of window to sample gaps from (window centered on current time)
        time_until_blocked: time until the channel is blocked (i.e. no more gaps)
        read_delay: delay between start of read and start of sequencing
    """
    def __init__(self, gaps_df, window_width, time_until_blocked, read_delay) -> None:
        # self.gaps_df = gaps_df
        assert set(gaps_df["gap_type"].cat.categories.values).issubset({"short", "long"})
        self.gaps_sorted_by_starts = gaps_df.sort_values("gap_start")
        self.gaps_sorted_by_ends = gaps_df.sort_values("gap_end")
        self.first_gap_start = gaps_df["gap_start"].min()
        self.last_gap_end = gaps_df["gap_end"].max()

        self.window_width = window_width
        self.time_until_blocked = time_until_blocked
        self.read_delay = read_delay
        
    def save(self, save_path):
        # gaps_sorted_by_starts enough
        np.savez_compressed(
            save_path, gaps_df=self.gaps_sorted_by_starts[["gap_start", "gap_end", "gap_duration", "gap_type"]].values, window_width=self.window_width, time_until_blocked=self.time_until_blocked, read_delay=self.read_delay,
        )
        
    @classmethod
    def load(cls, save_path):
        # loading pandas object
        data = np.load(save_path, allow_pickle=True)
        # workaround to load/save pd Dataframe with np.savez_compressed
        def modify_args(*, gaps_df, **kwargs):
            gaps_df = pd.DataFrame(gaps_df, columns=["gap_start", "gap_end", "gap_duration", "gap_type"])
            gaps_df["gap_type"] = gaps_df["gap_type"].astype("category")
            kwargs.update(gaps_df=gaps_df)
            return kwargs
        return cls(**modify_args(**data))
            
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None, long_gap_threshold=None, window_width=None):
        """mixes gaps from all channels"""
        if read_delay is None:
            read_delay = compute_median_read_delay(seqsum_df)
        
        if window_width is None:
            window_width = seqsum_df["end_time"].max() / 12
        
        dist = fit_seq_time_mle(seqsum_df)
        
        gaps_df = get_gaps_df_including_last_gap(seqsum_df)
        
        if long_gap_threshold is None:
            # single long gap threshold, not per window (for simplicity)
            long_gap_threshold = compute_long_gap_threshold(gaps_df["gap_duration"].values)
        gaps_df["gap_type"] = np.where(gaps_df["gap_duration"] <= long_gap_threshold, "short", "long")
        gaps_df["gap_type"] = gaps_df["gap_type"].astype("category")
        
        channel_nb = 0
        def get_gap_sampler(random_state=None):
            time_until_blocked = dist.rvs(random_state=random_state)
            nonlocal channel_nb
            channel_nb += 1
            return channel_nb, cls(gaps_df, window_width=window_width, time_until_blocked=time_until_blocked, read_delay=read_delay)
        
        return get_gap_sampler
        
    def __repr__(self):
        return f"RollingWindowGapSampler(window_width={self.window_width}, time_until_blocked={self.time_until_blocked}, read_delay={self.read_delay}, nb_gaps={len(self.gaps_sorted_by_starts)})"
    
    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        return self.read_delay
    
    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        t = channel_stats.time_active_without_mux_scan
        if (t >= self.time_until_blocked) or (len(self.gaps_sorted_by_starts) == 0):
            return GapSampler.GapType.Broken, np.inf
        
        # compute window to sample from
        window = np.array([t - self.window_width/2, t + self.window_width/2])
        # shift window if outside sequencing range so it has at least one gap
        if window[0] < self.first_gap_start:
            window = window + (self.first_gap_start - window[0])
        if window[1] > self.last_gap_end:
            window = window - (window[1] - self.last_gap_end)

        # find gaps in window, [i1, j1): all gaps starting in window, [i2, j2): all gaps ending in window, and gaps spanning the entire window
        # note: gaps that start and end in the window are twice as likely to be selected
        # could be made more efficient knowing that the time only increases, but log(n) binary search should be fast enough
        i1 = np.searchsorted(self.gaps_sorted_by_starts["gap_start"], window[0], side="left")
        j1 = np.searchsorted(self.gaps_sorted_by_starts["gap_start"], window[1], side="right")
        i2 = np.searchsorted(self.gaps_sorted_by_ends["gap_end"], window[0], side="left")
        j2 = np.searchsorted(self.gaps_sorted_by_ends["gap_end"], window[1], side="right")
        # i1, j1, i2, j2
        
        if j1-i1 + j2-i2 == 0:
            # gap spans entire window or no gap in window (long read)
            # assert i1 == j1 and i2 == j2
            assert i1 > 0
            gap_info = self.gaps_sorted_by_starts.iloc[i1-1]
            if gap_info["gap_start"] < window[0] and gap_info["gap_end"] > window[1]:
                # gap spans entire window
                return GapSampler.GapType.Long, gap_info["gap_duration"]
            else:
                # no gap in window (long read), so get surrounding gaps
            
                previous_gap_end = self.gaps_sorted_by_ends.iloc[i2-1]["gap_end"] if i2 > 0 else "na"
                next_gap_start = self.gaps_sorted_by_starts.iloc[j1]["gap_start"] if j1 < len(self.gaps_sorted_by_starts) else "na"
                # logger.info(f"no gaps in window {window}, {(i1, j1, i2, j2)}, previous gap end {previous_gap_end}, next gap start {next_gap_start}") # logger config not transferred to parallel process with joblib (i.e. not writing to stdout)
                print(f"no gaps starting or ending in window {window}, {(i1, j1, i2, j2)}, previous gap end {previous_gap_end}, next gap start {next_gap_start}")
            
                i1 = max(i1-10, 0)
                j1 = min(j1+10, len(self.gaps_sorted_by_starts))
                i2 = max(i2-10, 0)
                j2 = min(j2+10, len(self.gaps_sorted_by_ends))
            
            
        if random_state.random() < (j1-i1)/(j1-i1 + j2-i2):
            gap_duration, gap_type = self.gaps_sorted_by_starts.iloc[random_state.integers(i1, j1)][["gap_duration", "gap_type"]]
        else:
            gap_duration, gap_type = self.gaps_sorted_by_ends.iloc[random_state.integers(i2, j2)][["gap_duration", "gap_type"]]
        
        return (GapSampler.GapType.Short if gap_type == "short" else GapSampler.GapType.Long), gap_duration

# not inheriting from GapSampler because it's not a GapSampler itself
class RollingWindowGapSamplerPerChannel:
    """
    Similar to RollingWindowGapSampler, but samples the gaps for each channel separately, not mixing gaps between channels.
    
    For each channel, it samples the corresponding channel in the originak dataset
    
    Called rolling_window_per_channel in the paper
    
    """
    def __init__(self):
        pass
    
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None, long_gap_threshold=None, window_width=None, n_channels_full=None):
        """
        Gap sampler per channel, sampling gaps from the corresponding channel in the original dataset
        
        First sample the channel, then take the gaps, long_gap_threshold from this channel (as well as read_delay, window_width from this channel, if not provided).
        
        Args:
            n_channels_full: total number of channels including completely broken ones (cannot be inferred from the sequencing summary file because these channels contain no reads)
                this determines the probability of sampling a broken channel
        """
        
        gaps_df_all = get_gaps_df_including_last_gap(seqsum_df)
        channels = list(seqsum_df["channel"].unique())
        if n_channels_full is None:
            n_channels_full = len(channels)
        assert n_channels_full >= len(channels)
        
        # do it outside to avoid a SettingWithCopyWarning
        # per channel, single long gap threshold, not per window (for simplicity)
        def add_gap_type(gaps_df_channel):
            long_gap_threshold_channel = compute_long_gap_threshold(gaps_df_channel["gap_duration"].values) if long_gap_threshold is None else long_gap_threshold
            gaps_df_channel["gap_type"] = np.where(gaps_df_channel["gap_duration"] <= long_gap_threshold_channel, "short", "long")
            return gaps_df_channel
        gaps_df_all = gaps_df_all.groupby("channel", observed=True).apply(add_gap_type)
        gaps_df_all["gap_type"] = gaps_df_all["gap_type"].astype("category")
        
        channel_nb = 0
        def get_gap_sampler(random_state=np.random.default_rng(2)):
            nonlocal channel_nb
            channel_nb += 1
            
            channel_idx = random_state.integers(n_channels_full)
            if channel_idx >= len(channels):
                # broken channel
                return channel_nb, CompletelyBrokenGapSampler()
            channel = channels[channel_idx]
            
            # gaps_df, seqsum_df_channel are views, don't modify them since they modify the parent object, see Question3 in answer https://stackoverflow.com/a/53954986/3932263, or call .copy()
            gaps_df = gaps_df_all[gaps_df_all["channel"] == channel]
            seqsum_df_channel = seqsum_df[seqsum_df["channel"] == channel]
            
            read_delay_channel = compute_median_read_delay(seqsum_df_channel) if read_delay is None else read_delay
            
            seq_time = seqsum_df_channel["end_time"].max()
            window_width_channel = seq_time / 12 if window_width is None else window_width
            
            return channel_nb, RollingWindowGapSampler(gaps_df, window_width=window_width_channel, time_until_blocked=seq_time, read_delay=read_delay_channel)
        
        return get_gap_sampler