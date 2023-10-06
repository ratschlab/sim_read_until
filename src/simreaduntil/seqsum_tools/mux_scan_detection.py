"""
Detect and remove mux scans (for non-Flongles)
"""

from collections import Counter
from pathlib import Path
from typing import List, Tuple
import numpy as np
import pandas as pd
import tqdm
from simreaduntil.seqsum_tools.seqsum_preprocessing import sort_and_clean_seqsum_df

from simreaduntil.shared_utils.utils import is_sorted

def find_mux_scans(start_times, end_times, muxes, max_block_gap=1, max_intv_gap=20, min_mux_frac=0.95) -> List[Tuple[float, float]]:
    """
    Find mux scan windows when the channel is scanning all the wells (muxes) to determine the best well
    
    **Note**: This does not work for Flongle runs because they only have a single mux.
    
    The mux scans are characterized by the mux changing from 1, 2, 3, 4 in this order. 
    Although they should occur at 90 minutes intervals, this is not always exactly the case.
    We refer to a mux scan as a consecutive set of windows (of reads at most max_block_gap apart), where the dominant
    fraction of muxes (according to min_mux_frac) of window i is i.
    
    This method works best with many channels because muxes are performed simultaneously 
    across all channels and this way, we can detect whether all channels are scanning the same mux at a time.
    Assumes that the sequencing starts at time 0.
    
    This code was refactored from sim_utils.py in UNCALLED, mostly by documenting the code and renaming variables.
    The code does not seem very robust, but we don't want to modify it too much.
    
    Args:
        start_times: start times of reads
        end_times: end times of reads
        muxes: muxes, i.e. which pore in the channel is used
        max_block_gap: maximum distance between reads in windows that will be tested for mux scan segment (for a fixed mux 1, 2, 3, 4), 
            if two reads are more than max_block_gap apart, they are in different blocks
        
        max_intv_gap: if two mux scan windows are more than max_intv_gap apart, they are in different mux scan windows
        min_mux_frac: minimum fraction to classify a region as scanning a certain mux (1, 2, 3 or 4)
        
    Returns:
        2d np.array of shape (m, 2) of (start, end) times of mux scans
"""
    # pd.Series does not work because of indexing (probably)
    assert isinstance(start_times, np.ndarray)
    assert isinstance(end_times, np.ndarray)
    assert isinstance(muxes, np.ndarray)

    order = np.argsort(start_times)
    start_times = start_times[order] # start times of reads
    end_times = end_times[order] # end times of reads
    muxes = muxes[order] # muxes, i.e. which pore in the channel is used

    # first find blocks within which reads are at most max_block_gap apart
    
    active_blocks = list() # set of (t_start, t_end) describing active periods
    window_start = start_times[0] # start time of current window
    window_end = end_times[0] # end time of current window

    for (rst, ren) in zip(start_times[1:], end_times[1:]):
        if rst - window_end > max_block_gap:
            active_blocks.append( (window_start, window_end) )
            window_start = rst
            window_end = ren
        else:
            window_end = max(ren, window_end)
            #cur_window_end = ren
    active_blocks.append((window_start, window_end))

    # given the blocks, find mux scan windows
    # they are characterized by the mux changing from 1, 2, 3, 4 in this order, so we only record a mux scan window when there are exactly four elements with increasing muxes
    # for the mux=1 window, we merge it with the previous mux=1 window if they are less than max_intv_gap apart; if they are more apart, we discard the previous mux=1 window
    # finally, the mux window is the window preceding the first mux=1 window to the window following the end of the last mux=4 window
    # whenever two windows are more than max_intv_gap apart, they are in different mux scan windows
    
    mux_scan_windows = list() # list of mux scans, each entry is a scan of four elements
    scan_gaps = list() # list of (float, float) tuples denoting gaps to the previous and next window for the mux=1 window and mux=4 window
    cur_scan_segments = list() # running list of four (start, end) tuples, corresponding to muxes 1,2,3,4 respectively

    prev_window_end = 0
    for (window_start, window_end) in active_blocks:

        if len(cur_scan_segments) > 0 and window_start - cur_scan_segments[-1][1] > max_intv_gap:
            if len(cur_scan_segments) == 4:
                mux_scan_windows.append(cur_scan_segments)
                # todo2: bug? scan_gaps not appended to
                raise Exception("Should not happen")
            cur_scan_segments = list()

        mux_counts = Counter(muxes[(start_times >= window_start) & (start_times < window_end)]) # count muxes in current window across all channels and see if all channels are predominantly using one mux
        mux_counts = [(c, m) for (m, c) in mux_counts.items()]
        top_count, top_mux = max(mux_counts)

        if top_count / sum((c for (c, m) in mux_counts)) >= min_mux_frac:
            # if the window has a dominant mux (>= min_mux_frac = 0.95)

            if top_mux != 4 and len(cur_scan_segments) == 4:
                # out-of-order, but since we went from 1 to 4, we can finish the mux scan window
                mux_scan_windows.append(cur_scan_segments)

                scan_gaps.append((gap_before_scan,
                                  window_start - cur_scan_segments[-1][1]))

                cur_scan_segments = list()

            if len(cur_scan_segments) > 0 and top_mux == len(cur_scan_segments):
                # window has same mux as previous window
                if window_end - cur_scan_segments[-1][1] < max_intv_gap:
                    # merge windows by modifying scan window
                    cur_scan_segments[-1] = (cur_scan_segments[-1][0], window_end)

                elif top_mux == 1:
                    # if it is mux=1 and above max_intv_gap, discard the old window
                    cur_scan_segments[0] = (window_start, window_end)
                    gap_before_scan = window_start - prev_window_end # gap before mux=1 window
                # otherwise: don't take window into account, todo2: should discard window?

            elif top_mux-1 == len(cur_scan_segments):
                # window has mux one higher than previous window
                cur_scan_segments.append( (window_start, window_end) )
                if len(cur_scan_segments) == 1:
                    # mux=1 window in mux scan
                    gap_before_scan = window_start - prev_window_end

            else:
                # in a mux scan, the mux only increases sequentially
                # since this was not the case, discard current scan
                cur_scan_segments = list()
        else:
            if len(cur_scan_segments) == 4:
                mux_scan_windows.append(cur_scan_segments)
                scan_gaps.append((gap_before_scan,
                                  window_start - cur_scan_segments[-1][1]))
            cur_scan_segments = list()

        prev_window_end = window_end


    # enlargen windows to include gaps
    scans = list()
    for (segs, gaps) in zip(mux_scan_windows, scan_gaps):
        scans.append((segs[0][0]-gaps[0], segs[-1][1]+gaps[1]))

    return np.array(scans) if len(scans) > 0 else np.zeros(shape=(0, 2), dtype=np.float64)

def check_mux_scan_windows(seqsum_df, mux_scan_windows):
    """
    Check that each read is either entirely outside or entirely inside a mux scan, 
    """
    assert all(np.diff(mux_scan_windows.flatten()) > 0), "windows should have positive length"
    
    for (start, end) in tqdm.tqdm(mux_scan_windows, "Checking no read partially intersects mux scan window"):
        # check each element either is before mux scan, inside, or after
        assert all(    
            ((seqsum_df["start_time"] <= start) & (seqsum_df["end_time"] <= start)) |
            ((seqsum_df["start_time"] >= start) & (seqsum_df["end_time"] <= end)) |
            ((seqsum_df["start_time"] >= end) & (seqsum_df["end_time"] >= end))
        )
        
def plot_mux_scans(mux_scan_windows, ax):
    """
    Add mux scans as gray regions (axvspan)
    
    Args:
        mux_scan_windows: list of tuples (start_time, end_time)
        ax: matplotlib axis
    """
    for (start_time, end_time) in mux_scan_windows: # [(10, 20), (25,30)]
        ax.axvspan(start_time, end_time, facecolor='gray', alpha=0.5)
        ax.axvline(start_time, color="black")
        ax.axvline(end_time, color="black")

def _remove_reads_inside_mux_scans(seqsum_df, mux_scan_windows):
    """
    Remove reads that are entirely inside mux scans
    
    Args:
        seqsum_df: seqsum dataframe
        mux_scan_windows: list of tuples (start_time, end_time)
    """
    
    # build mask of elements that are not inside mux scans
    removal_mask = np.zeros(len(seqsum_df), dtype=bool)
    for (start, end) in mux_scan_windows:
        removal_mask |= (seqsum_df["start_time"] >= start) & (seqsum_df["end_time"] <= end)
    # remove reads inside a mux scan

    return seqsum_df[~removal_mask]
        
def remove_mux_scans(seqsum_df, mux_scan_windows):
    """
    Remove the mux scans from the sequencing summary by shifting times
    
    Removes reads that are entirely inside mux scans
    Adds a column "nb_scans_before", i.e. how many mux scans have been performed before this read
    
    Args:
        seqsum_df: sequencing summary
        mux_scan_windows: list of (start, end) tuples of mux scan windows
    
    Returns:
        Tuple of (the modified sequencing summary, mux scan boundaries)
    """
    
    seqsum_df = seqsum_df.copy()
    mux_scan_windows = mux_scan_windows.copy()
    
    seqsum_df = _remove_reads_inside_mux_scans(seqsum_df, mux_scan_windows)
    seqsum_df["nb_scans_before"] = 0
    
    assert is_sorted(mux_scan_windows.flatten()), "mux scan windows should be sorted"
    
    for i in range(len(mux_scan_windows)):
        start, end = mux_scan_windows[i]
        # remove ith mux scan from times
        window_length = end - start
        
        mux_scan_windows[mux_scan_windows > start] -= window_length
        
        seqsum_df.loc[seqsum_df["start_time"] > start, "start_time"] -= window_length
        seqsum_df.loc[seqsum_df["end_time"] > start, "end_time"] -= window_length
        seqsum_df.loc[seqsum_df["template_start"] > start, "template_start"] -= window_length
        
        seqsum_df.loc[seqsum_df["start_time"] > start, "nb_scans_before"] += 1
        
    assert np.allclose(mux_scan_windows[:, 0], mux_scan_windows[:, 1]), "mux scan windows should now have length 0"
    
    return seqsum_df, mux_scan_windows[:, 0]

def get_seqsum_filename_with_removed_mux_scans(seqsum_filename):
    """Filename to remove sequencing summary with removed mux scans, in same directory"""
    seqsum_filename = Path(seqsum_filename)
    return seqsum_filename.parent / ("no_mux_scans_" + seqsum_filename.name)

def remove_mux_scans_from_file(sequencing_summary_file, save_filename=None):
    """
    Remove mux scans from sequencing summary file and write it back to a new file
    
    Args:
        sequencing_summary_file: path to sequencing summary file
        save_filename: filename to save to; if None, don't save
        
    Returns:
        Tuple of (original sequencing summary, mux scan boundaries, path to new sequencing summary file)
    """
    
    df_read = pd.read_csv(sequencing_summary_file, sep="\t")
    seqsum_df = sort_and_clean_seqsum_df(df_read, min_gap_duration=None)#, min_gap_duration=0.05)

    mux_scan_windows = find_mux_scans(seqsum_df["start_time"].values, seqsum_df["end_time"].values, seqsum_df["mux"].values)
    check_mux_scan_windows(seqsum_df, mux_scan_windows)
    print("Removing mux scans")
    seqsum_df_new, mux_scan_boundaries = remove_mux_scans(seqsum_df, mux_scan_windows)
    
    if save_filename is not None:
        print(f"Saving to '{save_filename}'")
        seqsum_df_new.to_csv(save_filename, sep="\t", index=False)
    
    return seqsum_df, mux_scan_boundaries