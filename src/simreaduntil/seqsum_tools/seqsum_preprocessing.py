"""
Utilities to preprocess and clean the sequencing summary dataframe before parameter extraction
"""

import numpy as np
import pandas as pd
import scipy
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.utils import is_sorted


logger = setup_logger_simple(__name__)
"""module logger"""

def ensure_min_gaps_between_reads_single(df_single, min_gap_duration: int):
    """
    Ensure that the gap between reads is at least min_gap, shifting reads accordingly if this is not the case
    
    Real data sometimes has reads with negative gap (e.g. -1e-12 or -1e-4). This poses problems 
    with the simulator, so we correct for it.
    
    Args:
        df_single: dataframe for single channel, dataframe with columns "start_time", "duration", "channel", "template_start"
        min_gap_duration: minimum gap between reads, can be negative; to ensure gaps >= 0, choose 1e-12 instead of 0 as otherwise the gap can be negative due to floating point errors
        
    Returns:
        dataframe, the same if the gap is already present
    """
    assert is_sorted(df_single["start_time"].values), "df must be sorted by start_time"
    assert df_single["channel"].nunique() == 1
    assert all(df_single["template_start"] >= df_single["start_time"])
    assert df_single["start_time"].iloc[0] >= 0
    # assert all(df["template_start"] - df["start_time"] <= df["duration"])

    prev_end_times = np.concatenate(([0], (df_single["start_time"] + df_single["duration"]).iloc[:-1].values))
    shifts = np.maximum(prev_end_times + min_gap_duration - df_single["start_time"].values, 0)
    if len(shifts) > 0 and shifts.max() > 0:
        logger.info(f"Reads overlapping by {shifts.max()} (including min_gap {min_gap_duration})")
    else:
        return df_single
    cum_shifts = shifts.cumsum()
    
    df_single = df_single.copy()
    df_single["start_time"] += cum_shifts
    df_single["template_start"] += cum_shifts
    
    if "end_time" in df_single.columns:
        df_single["end_time"] = df_single["start_time"] + df_single["duration"]
    if "prev_gap_duration" in df_single.columns:
        df_single["prev_gap_duration"] += shifts
        
    return df_single

def ensure_min_gaps_between_reads(seqsum_df, min_gap_duration: int):
    """
    Wrapper around ensure_min_gaps_between_reads_single for multiple channels
    """
    return seqsum_df.groupby("channel", group_keys=False, observed=True).apply(lambda df: ensure_min_gaps_between_reads_single(df, min_gap_duration=min_gap_duration))

def sort_and_clean_seqsum_df(seqsum_df, min_gap_duration=0):
    """
    Preprocess the sequencing summary dataframe, add extra information
    """
    
    seqsum_df = seqsum_df.copy()
    
    # map channels to contiguous range [1, num_channels]
    channel_mapping = {channel: i + 1 for (i, channel) in enumerate(sorted(seqsum_df["channel"].unique()))}
    seqsum_df["channel"] = seqsum_df["channel"].map(channel_mapping)
    
    seqsum_df["end_time"] = seqsum_df["start_time"] + seqsum_df["duration"]
    seqsum_df.sort_values(by="end_time", inplace=True)
    
    if min_gap_duration is not None:
        seqsum_df = ensure_min_gaps_between_reads(seqsum_df, min_gap_duration=min_gap_duration)
    
    return seqsum_df

def add_previous_gap_duration(seqsum_df, seq_start_time=0, update_if_present=False):
    """
    Add gap durations before read
    
    It assumes that the dataframe is sorted by channel and time
    This only gets the gaps preceding a read, so the last gap after the last read 
    until the sequencing end is not included.
    
    Args:
        seqsum_df: dataframe with columns "start_time", "end_time", "channel"
        seq_start_time: start time of the whole sequencing run
        update_if_present: whether to update the column if it is already present
    """
    if "prev_gap_duration" in seqsum_df.columns and not update_if_present:
        # optimization because the groupby apply is slow
        return seqsum_df
    
    assert all(seqsum_df.groupby("channel", observed=True).apply(lambda df: is_sorted(df["start_time"]))), "dataframe is not sorted by channel and time"
    seqsum_df = seqsum_df.copy()
    
    # gap durations: duration of gap just before the read
    # use group_keys=False and index=df.index!, see https://gist.github.com/maximilianmordig/d396a58a2455087438ac25155020dbb5
    seqsum_df["prev_gap_duration"] = seqsum_df.groupby("channel", group_keys=False, observed=True).apply(
        lambda df: pd.Series(df["start_time"].values - np.concatenate(([seq_start_time], df["end_time"].values[:-1])), index=df.index)
    ).squeeze(axis="rows") # squeeze because pd groupby apply returns a dataframe instead of a series when the groupby column has only one unique value
    return seqsum_df

def get_gaps_single_channel(df_single, seq_start_time, check_single_channel=True, check_no_overlap=True):
    """
    Get (ordered) gaps, the channel always starts with a gap and all gaps until the last read are returned
    
    Args:
        df_single: dataframe with reads for a single channel (start_time and end_time columns)
        seq_start_time: start time of the whole sequencing run
        check_single_channel: whether to check that a single channel is present (requires channel column)
        check_no_overlap: whether to check that there is no overlap between reads
    Returns:
        (gap_durations, gap_starts)
    """
    if len(df_single) == 0:
        # return np.array([]), np.array([])
        raise ValueError("No reads for this channel")
    if check_single_channel:
        assert df_single["channel"].nunique() == 1
    assert is_sorted(df_single["start_time"])
    assert seq_start_time <= df_single["start_time"].min()

    if check_no_overlap:
        assert all(df_single["end_time"].iloc[:-1].values <= df_single["start_time"].iloc[1:].values), f"""reads overlap by up to {(df_single["end_time"].iloc[:-1].values - df_single["start_time"].iloc[1:].values).max()}, please correct this"""

    gap_starts = np.concatenate(([seq_start_time], df_single["end_time"].values[:-1]))
    return df_single["start_time"].values - gap_starts, gap_starts

def compute_long_gap_threshold(gap_durations):
    """
    Compute long_gap_threshold as quantile of all gaps (excluding last gap)
    
    Does not include the last gap because the channel is probably blocked at the end of a run
    """
    # long_gap_threshold = np.median(gap_durations) + np.std(gap_durations)
    
    # long_gap_threshold = np.quantile(gap_durations, q=0.95)
    # alternatives:
    # np.median(gap_durations) + 5*scipy.stats.median_abs_deviation(gap_durations)
    long_gap_threshold = np.median(gap_durations) + 5*scipy.stats.iqr(gap_durations)

    return long_gap_threshold

def compute_median_read_delay(seqsum_df):
    """
    Compute median read delay across all channels
    We observed that the delay is rather constant across channel and does not change too much across channels

    Args:
        seqsum_df: (preprocessed) dataframe from a sequencing summary file
    """
    return np.median(seqsum_df["template_start"] - seqsum_df["start_time"])

def compute_median_pore_speed(seqsum_df):
    """Compute median reading speed in bp/s"""
    pore_speed = seqsum_df["sequence_length_template"] / seqsum_df["template_duration"]
    return np.median(pore_speed)