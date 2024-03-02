"""
Sample short and long gaps in windows, and eventually block the channel. Gaps are extracted from a 
window surrounding the current window.
"""

from typing import Tuple
import numpy as np
import tqdm

from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_long_gap_threshold
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler, fit_nb_bps_until_blocked, fit_seq_time_beta, fit_seq_time_mle, get_gaps_df_including_last_gap, restrict_gaps_df_to_window
from simreaduntil.seqsum_tools.seqsum_preprocessing import compute_median_read_delay


def compute_time_and_aggregation_windows(seq_end_time, nb_windows, fraction_overlap=0.5):
    """
    Compute time windows of same length, and data aggregation windows overlapping by fraction_overlap on each side
    
    Args:
        seq_end_time: end time of the sequencing run, [0, seq_end_time]
        nb_windows: number of windows to divide the data into
        fraction_overlap: fraction of window_length overlapping on each side, so total overlap is 2*fraction_overlap*window_length;
            The first window is only overlapping by fraction_overlap on the right, the last window only on the left.
        
    Returns:
        time_windows: array of time windows (window_start, window_end) when windows are active
        aggregation_time_windows: array of time windows (window_start, window_end) when data is aggregated
    """
    assert nb_windows >= 1
    assert 0 <= fraction_overlap <= 1
    
    time_windows = np.array(np.lib.stride_tricks.sliding_window_view(np.linspace(0, seq_end_time, nb_windows + 1), window_shape=2))
    window_length = time_windows[0, 1] - time_windows[0, 0]
    overlap = fraction_overlap * window_length
    aggregation_time_windows = np.array([(max(0, start - overlap), min(end + overlap, seq_end_time)) for (start, end) in time_windows])
    return time_windows, aggregation_time_windows

def dict_without_items(d, keys_to_remove):
    # if d is a npzfile, it does not read d[k] if the key is in keys_to_remove
    return {k: d[k] for k in d if k not in keys_to_remove}
    
class GapSamplerPerWindowUntilBlocked(GapSampler):
    """
    Gap sampler that has separate sample distributions per time window and eventually blocks.
    
    Called window_all_channels in the paper
    
    Args:
        short_gaps_per_window: array of short gaps for each time window
        long_gaps_per_window: array of long gaps for each time window
        time_windows: array of time windows (window_start, window_end)
        time_until_blocked: number of bps until the channel is blocked (one it is at least this number)
        read_delay: delay between read starting and first bp being read
    """
    def __init__(self, short_gaps_per_window, long_gaps_per_window, time_windows, time_until_blocked, read_delay):
        self.short_gaps_per_window = short_gaps_per_window
        self.long_gaps_per_window = long_gaps_per_window
        self.time_windows = time_windows
        self.time_until_blocked = time_until_blocked
        self.read_delay = read_delay
        
        self.cur_window_idx = 0
        
    def save(self, save_path):
        # np cannot create a ragged array (i.e. array of subarrays with different lengths), so we concatenate them and save the offsets
        # alternative is to use pickle, which is slower, but more space-efficient
        short_gaps_per_window_concat = np.concatenate(self.short_gaps_per_window)
        short_gaps_per_window_offsets = np.cumsum([len(gaps) for gaps in self.short_gaps_per_window])
        long_gaps_per_window_concat = np.concatenate(self.long_gaps_per_window)
        long_gaps_per_window_offsets = np.cumsum([len(gaps) for gaps in self.long_gaps_per_window])
        np.savez_compressed(
            save_path, time_windows=self.time_windows, time_until_blocked=self.time_until_blocked, read_delay=self.read_delay,
            short_gaps_per_window_concat=short_gaps_per_window_concat, short_gaps_per_window_offsets=short_gaps_per_window_offsets,
            long_gaps_per_window_concat=long_gaps_per_window_concat, long_gaps_per_window_offsets=long_gaps_per_window_offsets
        )
        
    @classmethod
    def load(cls, save_path):
        data = np.load(save_path)
        short_gaps_per_window = np.split(data["short_gaps_per_window_concat"], data["short_gaps_per_window_offsets"])
        long_gaps_per_window = np.split(data["long_gaps_per_window_concat"], data["long_gaps_per_window_offsets"])
        
        # return cls(
        #     short_gaps_per_window=short_gaps_per_window, short_gaps_per_window=long_gaps_per_window, 
        #     time_windows=data["time_windows"], time_until_blocked=data["time_until_blocked"], read_delay=data["read_delay"]
        # )
        return cls(
            **dict_without_items(data, ["short_gaps_per_window_concat", "short_gaps_per_window_offsets", "long_gaps_per_window_concat", "long_gaps_per_window_offsets"]), 
            short_gaps_per_window=short_gaps_per_window, long_gaps_per_window=long_gaps_per_window
        )
    
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None, time_and_aggregation_windows=None):
        """
        Compute parameters from a sequencing summary dataframe
        
        Call add_previous_gap_duration before calling this function
        It assumes that each channel becomes broken after its final read.
        
        Args:
            seqsum_df: sequencing summary dataframe
            read_delay: delay between read starting and first bp being read; if None, compute median read delay
            time_and_aggregation_windows: tuple of array of time windows (window_start, window_end) and array 
                of data aggregation windows (window_start, window_end);
                if None, use 4h windows with 50% overlap, i.e. [t, t+4] window with data from [t-2, t+6]
        
        Returns:
            function to create a gap sampler, so it is flexible with respect to the number of channels
        """
        if read_delay is None:
            read_delay = compute_median_read_delay(seqsum_df)
        
        # dist = fit_nb_bps_until_blocked(seqsum_df) # leads to a bimodal distribution of sequencing time per channel (histogram), but the real one is unimodal
        # dist = fit_seq_time_beta(seqsum_df)
        dist = fit_seq_time_mle(seqsum_df)

        if time_and_aggregation_windows is None:
            seq_end_time = seqsum_df["end_time"].max()
            nb_windows = max(1, round(seq_end_time / (4 * 3600))) # 1 window roughly every 4h
            time_windows, aggregation_time_windows = compute_time_and_aggregation_windows(seq_end_time, nb_windows=nb_windows, fraction_overlap=0.5)
        else:
            time_windows, aggregation_time_windows = time_and_aggregation_windows

        gaps_df = get_gaps_df_including_last_gap(seqsum_df)
        short_gaps_per_window = []
        long_gaps_per_window = []
        for (window_start, window_end) in tqdm.tqdm(aggregation_time_windows, desc="Computing gaps per window"):
            restricted_gaps = restrict_gaps_df_to_window(gaps_df, window_start, window_end)["gap_duration"].values
            long_gap_threshold = compute_long_gap_threshold(restricted_gaps) # recompute long gaps threshold!
            short_gaps_per_window.append(restricted_gaps[restricted_gaps <= long_gap_threshold])
            long_gaps_per_window.append(restricted_gaps[restricted_gaps > long_gap_threshold])
        
        channel_nb = 0
        def get_gap_sampler(random_state=None):
            time_until_blocked = dist.rvs(random_state=random_state)
            nonlocal channel_nb
            channel_nb += 1
            return channel_nb, cls(short_gaps_per_window, long_gaps_per_window, time_windows, time_until_blocked=time_until_blocked, read_delay=read_delay)
        
        return get_gap_sampler
    
    def __repr__(self):
        return f"GapSamplerPerWindowUntilBlocked(nb short_gaps_per_window={len(self.short_gaps_per_window)}, nb long_gaps_per_window={len(self.long_gaps_per_window)}, time_windows={self.time_windows}, time_until_blocked={self.time_until_blocked}, read_delay={self.read_delay})"
    
    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        return self.read_delay
    
    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        t = channel_stats.time_active_without_mux_scan
        while (self.cur_window_idx < len(self.time_windows)) and (t > self.time_windows[self.cur_window_idx][1]):
            self.cur_window_idx += 1
        # if (channel_stats.reads.number_bps_read >= self.bps_until_blocked) or (self.cur_window_idx >= len(self.time_windows)):
        if (t >= self.time_until_blocked) or (self.cur_window_idx >= len(self.time_windows)):
            return GapSampler.GapType.Broken, np.inf
        
        nb_short_gaps = len(self.short_gaps_per_window[self.cur_window_idx])
        nb_long_gaps = len(self.long_gaps_per_window[self.cur_window_idx])
        
        if isinstance(random_state, np.random.Generator):
            gap_idx = random_state.integers(nb_short_gaps + nb_long_gaps)
        else:
            gap_idx = random_state.randint(nb_short_gaps + nb_long_gaps)
            
        if gap_idx < nb_short_gaps:
            return GapSampler.GapType.Short, self.short_gaps_per_window[self.cur_window_idx][gap_idx]
        else:
            return GapSampler.GapType.Long, self.long_gaps_per_window[self.cur_window_idx][gap_idx - nb_short_gaps]
    