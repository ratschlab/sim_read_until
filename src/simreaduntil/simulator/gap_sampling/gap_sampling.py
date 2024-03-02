"""
Sample gaps and read delays in the channel.
"""

import enum
from typing import Any, Dict, Generator, Tuple
import numpy as np
import pandas as pd
import scipy

from simreaduntil.shared_utils.dna import get_random_DNA_seq
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.simulator.channel_stats import ChannelStats

logger = setup_logger_simple(__name__)
"""module logger"""

class GapSampler:
    """
    Sample gaps for a single channel
    """
    
    def set(*args, **kwargs):
        """
        Set attributes on this channel
        """
        raise NotImplementedError()
    
    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        """
        Delay between a read starting and the first bp being read. 
        
        Negative values returned by this function are treated as zero.
        
        Args:
            channel_stats: channel state
            random_state: random state to use for sampling
        """
        raise NotImplementedError()
    
    class GapType(str, enum.Enum):
        Short = "short_gap" # pore blocked for short amount of time, ended when a mux scan starts
        Long = "long_gap" # pore blocked for longer amount of time, interrupted during mux scan
        Broken = "broken_channel" # broken channel
        
    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapType, float]:
        """
        Time between reads.
        
        Args:
            channel_stats: channel state
            random_state: random state to use for sampling
        
        Returns:
            (gap_type, gap_duration)
        """
        raise NotImplementedError()
    
    def mark_long_gap_end(self, channel_stats: ChannelStats):
        """
        Notify about when the long gap ends. 
        
        This is useful when mux scans interrupt the long gap to avoid taking the time for them into account.
        
        Args:
            channel_stats: channel state
        """
        pass
    
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None):
        """
        Extract parameters for this class from a sequencing summary dataframe
        
        Args:
            seqsum_df: sequencing summary dataframe
            read_delay: delay between read starting and first bp being read
            
        Returns:
            channel_nb consecutive starting at 1, gap_sampler as an instance of the class cls
        """
        raise NotImplementedError()
    
    def save(self, save_path):
        """
        Save gap sampler to a file
        
        It should be saved so that the np.load function works with the constructor.
        
        Args:
            save_path: where to save
        """
        raise NotImplementedError()
    
    @classmethod
    def load(cls, save_path):
        """
        Load gap sampler from a file
        
        Args:
            save_path: where to load from
        """
        return cls(**np.load(save_path, allow_pickle=True))
    
class RandomGapSampler(GapSampler):
    """
    Gap sampler with random gap lengths, for testing mostly
    
    Channel breaks with some probability
    
    Called random_gaps in the paper
    """
    def __init__(self, prob_long_gap=0.5) -> None:
        super().__init__()
        
        self.prob_long_gap = prob_long_gap
        
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None):
        """Completely ignores the provided arguments"""
        channel_nb = 0
        def wrapper(random_state=None):
            nonlocal channel_nb
            channel_nb += 1
            return channel_nb, cls()
        return wrapper
    
    def save(self, save_path):
        np.savez_compressed(save_path, prob_long_gap=self.prob_long_gap)
        
    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        return random_state.uniform(0.1, 0.2)
    
    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        r = random_state.uniform()
        if r < self.prob_long_gap:
            return (self.GapType.Long, random_state.uniform(7, 10))
        elif r < 0.99:
            return (self.GapType.Short, random_state.uniform(0.1, 0.4))
        else:
            return (self.GapType.Broken, 0)
        
class CompletelyBrokenGapSampler(GapSampler):
    """
    Gap sampler for a broken channel, i.e. immediately broken
    """
    @classmethod
    def from_seqsum_df(cls, seqsum_df, read_delay=None):
        """Completely ignores the provided arguments"""
        channel_nb = 0
        def wrapper(random_state=None):
            nonlocal channel_nb
            channel_nb += 1
            return channel_nb, cls()
        return wrapper
    
    def save(self, save_path):
        np.savez_compressed(save_path, class_name=str(type(self)))
    
    @classmethod
    def load(cls, save_path):
        data = np.load(save_path)
        assert data["class_name"] == str(cls)
        return cls()
        
    def sample_read_start_delay(self, channel_stats: ChannelStats, random_state) -> float:
        logger.warning("Should never be called because channel is broken")
        return 0
    
    def sample_next_gap(self, channel_stats: ChannelStats, random_state) -> Tuple[GapSampler.GapType, float]:
        return (self.GapType.Broken, 0)
        
def compute_prob_long_gap(time_short_gap, time_long_gap, time_in_long_gaps_fraction):
    """
    Compute the probability to pick long gap to ensure that a fraction of time is spent in long gaps.
    
    The time_short_gap and time_long_gap can be the expected values. Then, this method returns
    the probability to draw long gaps such that the expected time spent in long gaps is 
    time_in_long_gaps_fraction of the total expected time.
    
    We take the median of short and long gaps, so the distribution is not preserved. This function 
    ensures that the fraction of time spent in long gaps is (approximately) preserved.
    
    Since time_long_gap >= time_short_gap, this decreases the probability of choosing a long gap.
    
    Args:
        time_short_gap: time spent in short gaps
        time_long_gap: time spent in long gaps
        time_in_long_gaps_fraction: fraction of time spent in long gaps
    """
    assert time_long_gap >= time_short_gap
    return 1 / (1 + (time_long_gap / time_short_gap) * (1 / time_in_long_gaps_fraction - 1))

def compute_shortgap_longgap_longgapprob(gap_durations, long_gap_threshold):
    """
    Compute the median of short and long gaps and the fraction of long gaps.
    """
    if len(gap_durations) == 0:
        return 0, np.inf, 1.0
    short_gap_median = np.median(gap_durations[gap_durations <= long_gap_threshold])
    long_gap_median = np.median(gap_durations[gap_durations > long_gap_threshold])
    # long_gap_fraction = sum(gap_durations > long_gap_threshold) / seqsum_df.shape[0]
    time_in_long_gaps_fraction = sum(gap_durations[gap_durations > long_gap_threshold]) / sum(gap_durations)
    return short_gap_median, long_gap_median, compute_prob_long_gap(short_gap_median, long_gap_median, time_in_long_gaps_fraction)

def fit_nb_bps_until_blocked(seqsum_df):
    """
    Fit a distribution to the number of bps until the channel is blocked
    """
    nb_bps_per_channel = seqsum_df.groupby("channel", observed=True)["sequence_length_template"].sum()
    dist = scipy.stats.weibull_min
    params = fit_distribution_and_check(dist, vals=nb_bps_per_channel)
    
    return dist(*params)

def fit_seq_time_beta(seqsum_df):
    """
    Fit a distribution to the sequencing end time of the channels
    
    Returns:
        distribution object
    """
    end_times = seqsum_df.groupby("channel", observed=True)["end_time"].max()
    dist = scipy.stats.beta
    params = fit_distribution_and_check(dist, vals=end_times)
    
    return dist(*params)

def fit_seq_time_mle(seqsum_df):
    """
    Fit MLE to sequencing times (i.e. discrete distribution that samples from it)
    
    Returns:
        distribution object
    """
    end_times = seqsum_df.groupby("channel", observed=True)["end_time"].max()
    return scipy.stats.rv_discrete(name="sequencing_time_distribution", values=(end_times, np.ones_like(end_times)/len(end_times)))

def fit_distribution_and_check(dist, vals):
    """
    Fit distribution to vals and check that p-value is >0.05
    
    Returns:
        fitted params
    """
    params = dist.fit(vals)
    
    # xvals = np.linspace(min(end_times), max(end_times), 100)
    # plt.plot(xvals, dist.pdf(xvals, *params))
    # plt.hist(end_times, bins=100, density=True) # normalized to sum to 1, not the same as pdf
    
    logger.info("Computing goodness of fit")
    try:
        goodness_res = scipy.stats.goodness_of_fit(dist, vals, n_mc_samples=100) # preferably use 9999, but takes more time
        logger.info(f"Goodness of fit result (p-value should be >0.05): {goodness_res}")
        if goodness_res.pvalue < 0.05:
            # null hypothesis: sample comes from the family of distributions
            # p-value < 0.05 -> null hypothesis should be rejected
            logger.warning(f"Goodness of fit test failed with p-value {goodness_res.pvalue}")
    except scipy.stats.FitError:
        # goodness-of-fit failed, probably due to repeated data due to too little data?
        logger.exception("Error computing goodness-of-fit")
        
    return params

def restrict_gaps_df_to_window(gaps_df, window_start, window_end):
    """
    Filter to gaps that are intersecting with the window [window_start, window_end]
    """
    mask = (
        ((window_start <= gaps_df["gap_start"]) & (gaps_df["gap_start"] <= window_end)) # gap start in window
        | ((window_start <= gaps_df["gap_end"]) & (gaps_df["gap_end"] <= window_end)) # gap end in window
        | ((gaps_df["gap_start"] <= window_start) & (window_end <= gaps_df["gap_end"])) # gap start before and gap end after window
    )
    return gaps_df[mask]

def get_gaps_df_including_last_gap(seqsum_df):
    """
    Get gaps in a dataframe, including last gap until sequencing end
    """
    gaps_df = pd.DataFrame({"gap_duration": seqsum_df["prev_gap_duration"], "gap_end": seqsum_df["start_time"], "channel": seqsum_df["channel"]})
    
    # add gap after last read, per channel
    seq_end_time = seqsum_df["end_time"].max()
    last_gaps_df = seqsum_df.groupby("channel", observed=True)["end_time"].max().to_frame(name="gap_start").reset_index()
    last_gaps_df["gap_end"] = seq_end_time
    last_gaps_df["gap_duration"] = seq_end_time - last_gaps_df["gap_start"]
    last_gaps_df = last_gaps_df[last_gaps_df["gap_duration"] > 0]
    
    gaps_df = pd.concat((gaps_df, last_gaps_df)).reset_index(drop=True)
    gaps_df["gap_start"] = gaps_df["gap_end"] - gaps_df["gap_duration"]
    
    return gaps_df
