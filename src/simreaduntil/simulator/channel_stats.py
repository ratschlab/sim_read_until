"""
Tracking channel statistics, for a single or multiple channels
"""

import threading
from typing import Iterable, List
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd
from simreaduntil.shared_utils.plotting import FIGURE_EXT
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.merge_axes import save_fig_and_pickle
from simreaduntil.shared_utils.plotting import make_tight_layout
from simreaduntil.shared_utils.plotting import _disable_x_ticks

from simreaduntil.simulator.utils import format_percentage

__all__ = ["ChannelStats", "combine_channel_stats"]

logger = setup_logger_simple(__name__)
"""module logger"""

class ElementStats:
    """
    Per element type, keep track of time spent, current number, finished number of elements
    
    Args:
        cur_number: current number of active elements
        finished_number: number of elements that have finished
        time_spent: time spent on elements (including active elements)
    """
    def __init__(self, cur_number=0, finished_number=0, time_spent=0):
        self.cur_number = cur_number
        self.finished_number = finished_number
        self.time_spent = time_spent
        
    def __repr__(self):
        return f"{self.__class__.__name__}(cur_number={self.cur_number}, finished_number={self.finished_number}, time_spent={self.time_spent})"
        
    def __eq__(self, other: "ElementStats"):
        return self.cur_number == other.cur_number and self.finished_number == other.finished_number and np.allclose(self.time_spent, other.time_spent)
    
    def check_consistent(self, n_channels_running):
        """Check whether the stats are consistent"""
        return 0 <= self.cur_number <= n_channels_running and 0 <= self.finished_number and 0 <= self.time_spent
    
    def combine_with(self, other: "ElementStats"):
        """
        Combine with another ElementStats object
        """
        self.cur_number += other.cur_number
        self.finished_number += other.finished_number
        self.time_spent += other.time_spent
        
    def start(self):
        """
        Start a new element
        """
        self.cur_number += 1
        
    def add_time(self, time, **kwargs):
        """
        Add time to an active element
        
        time can be negative
        """
        self.time_spent += time
        
    def start_and_add_time(self, time, **kwargs):
        """Combines start and add_time"""
        self.start()
        self.add_time(time, **kwargs)
        
    def finish(self, **kwargs):
        """Finish the active element"""
        self.finished_number += 1
        self.cur_number -= 1
        
    def add_time_and_finish(self, time, **kwargs):
        """Combines start and add_time"""
        self.add_time(time, **kwargs)
        self.finish(**kwargs)
        
    def add_full(self, time, **kwargs):
        """Start, add_time and finish"""
        self.start()
        self.add_time(time, **kwargs)
        self.finish(**kwargs)
        
class ReadElementStats(ElementStats):
    """
    Extends ElementStats with read-specific statistics
    """
    def __init__(self, cur_number=0, finished_number=0, time_spent=0, 
                 number_bps_requested=0, number_bps_read=0, number_bps_rejected=0,
                 cur_number_stop_receiving=0, fin_number_stop_receiving=0, 
                 fin_number_rejected=0,
                 number_rejected_missed=0, number_stop_receiving_missed=0):
        super().__init__(cur_number=cur_number, finished_number=finished_number, time_spent=time_spent)
        
        self.cur_number_stop_receiving = cur_number_stop_receiving
        self.fin_number_stop_receiving = fin_number_stop_receiving # number of reads that were set to stop_receiving when they finished
        self.fin_number_rejected = fin_number_rejected # number of reads that were successfully rejected, includes mux scan unblock etc
        
        self.number_bps_requested = number_bps_requested # number of basepairs that were requested via ReadUntil
        self.number_bps_read = number_bps_read # number of basepairs that were read
        self.number_bps_rejected = number_bps_rejected # number of basepairs rejected by ReadUntil
        
        # missed actions
        self.number_rejected_missed = number_rejected_missed # number of reads that were rejected, but that were no longer current when the decision was received
        self.number_stop_receiving_missed = number_stop_receiving_missed # number of reads that were set to stop receiving, but that were no longer current when the decision was received or were rejected afterwards
        
    def __repr__(self):
        return f"ReadElementStats(cur_number={self.cur_number}, finished_number={self.finished_number}, time_spent={self.time_spent}, number_bps_requested={self.number_bps_requested}, number_bps_read={self.number_bps_read}, number_bps_rejected={self.number_bps_rejected}, cur_number_stop_receiving={self.cur_number_stop_receiving}, fin_number_stop_receiving={self.fin_number_stop_receiving}, fin_number_rejected={self.fin_number_rejected}, number_rejected_missed={self.number_rejected_missed}, number_stop_receiving_missed={self.number_stop_receiving_missed})"
    
    def __eq__(self, other: "ReadElementStats"):
        return super().__eq__(other) and \
            self.cur_number_stop_receiving == other.cur_number_stop_receiving and self.fin_number_stop_receiving == other.fin_number_stop_receiving \
                and self.fin_number_rejected == other.fin_number_rejected and self.number_bps_requested == other.number_bps_requested \
                    and self.number_bps_read == other.number_bps_read and self.number_bps_rejected == other.number_bps_rejected \
                        and self.number_rejected_missed == other.number_rejected_missed and self.number_stop_receiving_missed == other.number_stop_receiving_missed
                        
    def check_consistent(self, n_channels_running):
        super().check_consistent(n_channels_running)
        assert 0 <= self.cur_number_stop_receiving <= self.cur_number
        assert 0 <= self.fin_number_stop_receiving <= self.finished_number
        assert 0 <= self.fin_number_rejected <= self.finished_number
        
        assert 0 <= self.number_bps_read
        assert 0 <= self.number_bps_rejected
        assert 0 <= self.number_bps_requested <= self.number_bps_read
        
        assert 0 <= self.number_rejected_missed
        assert 0 <= self.number_stop_receiving_missed
        
    def combine_with(self, other: "ReadElementStats"):
        super().combine_with(other)
        self.cur_number_stop_receiving += other.cur_number_stop_receiving
        self.fin_number_stop_receiving += other.fin_number_stop_receiving
        self.fin_number_rejected += other.fin_number_rejected
        self.number_bps_read += other.number_bps_read
        self.number_bps_rejected += other.number_bps_rejected
        
        self.number_bps_requested += other.number_bps_requested
        self.number_rejected_missed += other.number_rejected_missed
        self.number_stop_receiving_missed += other.number_stop_receiving_missed
        
    def add_time(self, time, **kwargs):
        super().add_time(time)
        if "nb_bps" in kwargs:
            self.number_bps_read = kwargs["nb_bps"]
        else:
            self.number_bps_read += kwargs["nb_new_bps"]
        
    def finish(self, **kwargs):
        super().finish(**kwargs)
        if kwargs["stopped_receiving"]:
            self.cur_number_stop_receiving -= 1
            self.fin_number_stop_receiving += 1
        nb_bps_rejected = kwargs.get("nb_bps_rejected", 0)
        if nb_bps_rejected > 0:
            self.number_bps_rejected += nb_bps_rejected
            self.fin_number_rejected += 1

class ChannelStats:
    """
    Keeps track of channel statistics for a single / multiple channels
    
    Tracked stats include:
        - time, finished number and current number of short and long gaps, mux scans, unblock delays, reads
        - for reads:     
            - number of successful and unsuccessful (if decision was too late) read actions (unblock, stop receiving)
            - number of bases read (not necessarily all returned), bases rejected, bases stopped receiving (but still read)

    Some statistics only concern finished elements (fin_*), others concern elements in-progress (cur_* and number_bps_read)
    
    Args:
        n_channels: number of channels, not necessarily all running
        short_gaps: ElementStats for short gaps
        long_gaps: ElementStats for long gaps
        unblock_delays: ElementStats for unblock delays
        mux_scans: ElementStats for mux scans
        channel_broken: ElementStats for channel broken
        reads: ReadElementStats for reads
        no_reads_left: ElementStats for no reads left
    """
    def __init__(self, 
                 n_channels: int=0,
                 short_gaps: ElementStats=None, long_gaps: ElementStats=None, unblock_delays: ElementStats=None, mux_scans: ElementStats=None,
                 channel_broken: ElementStats=None, reads: ReadElementStats=None, no_reads_left: ElementStats=None
                 ):
        self.n_channels = n_channels # number of channels, not necessarily all running
          
        self.short_gaps = ElementStats() if short_gaps is None else short_gaps
        self.long_gaps = ElementStats() if long_gaps is None else long_gaps
        self.unblock_delays = ElementStats() if unblock_delays is None else unblock_delays
        self.mux_scans = ElementStats() if mux_scans is None else mux_scans
        self.reads = ReadElementStats() if reads is None else reads
        self.no_reads_left = ElementStats() if no_reads_left is None else no_reads_left
        self.channel_broken = ElementStats() if channel_broken is None else channel_broken
        
        self.check_consistent()
      
    def combine_with(self, other: "ChannelStats"):
        """
        Combine with another ChannelStats object, e.g. for combining channel stats
        """
        self.n_channels += other.n_channels
                
        self.short_gaps.combine_with(other.short_gaps)
        self.long_gaps.combine_with(other.long_gaps)
        self.unblock_delays.combine_with(other.unblock_delays)
        self.mux_scans.combine_with(other.mux_scans)
        self.reads.combine_with(other.reads)
        self.no_reads_left.combine_with(other.no_reads_left)
        self.channel_broken.combine_with(other.channel_broken)
        
        self.check_consistent()
    
    def __eq__(self, other: "ChannelStats"):
        return self.n_channels == other.n_channels \
            and self.short_gaps == other.short_gaps and self.long_gaps == other.long_gaps \
                and self.unblock_delays == other.unblock_delays and self.mux_scans == other.mux_scans \
                    and self.reads == other.reads and self.no_reads_left == other.no_reads_left \
                        and self.channel_broken == other.channel_broken
                    
    def check_consistent(self, n_channels_running=None):
        """
        Checks whether statistics stored in this object make sense
        
        n_channels_running: to overwrite, e.g. if channel is stopped
        """
        if n_channels_running is None:
            n_channels_running = self.n_channels_running
                    
        self.short_gaps.check_consistent(n_channels_running)
        self.long_gaps.check_consistent(n_channels_running)
        self.unblock_delays.check_consistent(n_channels_running)
        self.mux_scans.check_consistent(n_channels_running)
        self.reads.check_consistent(n_channels_running)
        self.no_reads_left.check_consistent(n_channels_running)
        self.channel_broken.check_consistent(n_channels_running)
        
        # not equal since a mux scan adds a long gap twice (by splitting it) and simulation may end right after gap
        assert self.reads.finished_number + self.reads.cur_number <= self.short_gaps.finished_number + self.long_gaps.finished_number, "a gap always comes before a read"
        assert self.long_gaps.finished_number >= 0
        assert self.unblock_delays.finished_number + self.unblock_delays.cur_number <= self.reads.fin_number_rejected, "when a read is successfully rejected, an unblock delay follows which is either done or in-progress" # <= since it may be stopped right at this moment
        
    def number_full_reads(self):
        """Number of fully sequenced reads (without reject and reads in-progress with stop_receiving)"""
        return self.reads.finished_number - self.reads.fin_number_rejected
    
    @property
    def n_channels_running(self):
        """Number of channels currently running (i.e. not stopped), including broken channels"""
        return self.reads.cur_number + self.short_gaps.cur_number + self.long_gaps.cur_number + self.unblock_delays.cur_number + \
            self.no_reads_left.cur_number + self.mux_scans.cur_number + self.channel_broken.cur_number
            
    @property
    def time_active(self):
        """Time that channel is active until current timepoint, excluding time without reads and when channel was broken"""
        return self.short_gaps.time_spent + self.long_gaps.time_spent + self.unblock_delays.time_spent + self.mux_scans.time_spent + self.reads.time_spent
    
    @property
    def time_active_without_mux_scan(self):
        return self.time_active - self.mux_scans.time_spent
    
    def fraction_reading_time(self):
        """Fraction of time that a channel is reading until current timepoint"""
        return self.reads.time_spent / self.time_active
    
    def fraction_mux_scan(self):
        """Fraction of time spent with mux scan until current timepoint"""
        return self.mux_scans.time_spent / self.time_active
    
    def fraction_short_gaps(self):
        """Fraction of time of short gaps (not counting long gaps) until current timepoint"""
        return self.short_gaps.time_spent / self.time_active
    
    def fraction_unblocking(self):
        """Fraction of time spent due to read unblocking until current timepoint"""
        return self.unblock_delays.time_spent / self.time_active
    
    def fraction_long_gaps(self):
        """Fraction of time of long gaps (not counting short gaps) until current timepoint"""
        return self.long_gaps.time_spent / self.time_active
    
    def __str__(self):
        width = 7 # display width for each number
        if self.time_active == 0:
            # was not running, avoid division by zero
            return f"Summed running time for {self.n_channels} ({self.n_channels_running} active) channel(s): {self.time_active:.2f} seconds; no time elapsed since channel started"
        
        return (
            f"Summed simulation time for {self.n_channels_running} channel(s): {self.time_active:.2f} seconds, {format_percentage(self.fraction_reading_time())}% reading time"
            f", {format_percentage(self.fraction_short_gaps())}% short gap time, {format_percentage(self.fraction_unblocking())}% unblocking time, {format_percentage(self.fraction_mux_scan())}% mux scan time"
            f", {format_percentage(self.fraction_long_gaps())}% long gap time, no reads for {self.no_reads_left.time_spent:.2f} seconds"
            "\n"
            f"{self.reads.finished_number:>{width}} reads (+ currently {self.reads.cur_number:>{width}}): "
            f"{self.number_full_reads():>{width}} fully read; "
            f"{self.reads.fin_number_rejected:>{width}} rejected ({self.reads.number_rejected_missed:>{width}} missed); "
            f"{self.reads.fin_number_stop_receiving:>{width}} stopped to receive ({self.reads.number_stop_receiving_missed:>{width}} missed, + currently {self.reads.cur_number_stop_receiving:>{width}}); "
            "\n"
            f"{self.reads.number_bps_read:>{width}} bps read; {self.reads.number_bps_requested:>{width}} bps processed; {self.reads.number_bps_rejected:>{width}} bps rejected; "
            "\n"
            f"{self.short_gaps.finished_number:>{width}} ShortGaps (+ currently {self.short_gaps.cur_number:>{width}}); "
            f"{self.long_gaps.finished_number:>{width}} LongGaps (+ currently {self.long_gaps.cur_number:>{width}}); "
            f"{self.unblock_delays.finished_number:>{width}} UnblockDelayGaps (+ currently {self.unblock_delays.cur_number:>{width}})"
            f"{self.mux_scans.finished_number:>{width}} MuxScans (+ currently {self.mux_scans.cur_number:>{width}})"
        )
        
    @staticmethod
    def get_table_header() -> str:
        """Get tab-separated header line for table output, see get_table_line"""
        header_names = [
                "active_time", "noreads_time", 
                "time_short_gaps", "time_long_gaps", "time_unblock_delays", "time_mux_scan",
                "num_reads",         "num_stop_receiving",     "num_short_gaps",     "num_long_gaps",     "num_unblock_delays", "num_rejected", "num_missed_rejected", "num_missed_stop_receiving",
                "num_mux_scans",
                "num_bases_read", "num_bases_rejected", "num_bases_processed",
                "cur_num_reads", "cur_num_stop_receiving", "cur_num_short_gaps", "cur_num_long_gaps", "cur_num_unblock_delays",
                "cur_num_mux_scans"
                "n_active_channels",
        ]
        return "\t".join(header_names)
    
    def get_table_line(self, header=False) -> str:
        """Get tab-separated line of stats formatted to be output into a file for plotting later on"""
        # keep in sync with get_table_header
        line = [
            self.time_active, self.no_reads_left.time_spent,
            self.short_gaps.time_spent, self.long_gaps.time_spent, self.unblock_delays.time_spent, self.mux_scans.time_spent,
            self.reads.finished_number, self.reads.reads.fin_number_stop_receiving, self.short_gaps.finished_number, self.long_gaps.finished_number, self.unblock_delays.finished_number, self.reads.fin_number_rejected, self.reads.number_rejected_missed, self.reads.number_stop_receiving_missed,
            self.mux_scans.finished_number,
            self.reads.number_bps_read, self.reads.number_bps_rejected, self.reads.number_bps_requested,
            self.reads.cur_number, self.reads.cur_number_stop_receiving, self.short_gaps.cur_number, self.long_gaps.cur_number, self.unblock_delays.cur_number,
            self.mux_scans.cur_number,
            self.n_channels_running,
        ]
        # width = 30
        # " ".join([f"{val:>{width}}" for val in line])
        return "\t".join(map(str, line))
    
def combine_stats(stats: Iterable[ChannelStats]):
    """aggregate stats from multiple channels into one"""
    res = ChannelStats()
    for stat in stats:
        res.combine_with(stat)
    return res

####################
# plotting
####################

def channel_stats_to_df(channel_stats: List[ChannelStats]):
    """
    Convert channel stats into a data.frame in wide format with one line per channel, useful for plotting
    """
    
    # check cur_number is zero for all elements
    assert all(all([
        chan.short_gaps.cur_number == 0, chan.long_gaps.cur_number == 0, chan.unblock_delays.cur_number == 0, chan.mux_scans.cur_number == 0, 
        chan.channel_broken.cur_number == 0, chan.reads.cur_number == 0, chan.no_reads_left.cur_number == 0
    ]) for chan in channel_stats), "cur number not zero"
    assert all(chan.n_channels == 1 for chan in channel_stats), "each element should represent one channel"

    df = pd.DataFrame(
        [(
            channel,
            channel.short_gaps.finished_number, channel.short_gaps.time_spent, 
            channel.long_gaps.finished_number, channel.long_gaps.time_spent, 
            channel.unblock_delays.finished_number, channel.unblock_delays.time_spent, 
            channel.mux_scans.finished_number, channel.mux_scans.time_spent, 
            channel.channel_broken.finished_number, channel.channel_broken.time_spent, 
            channel.reads.finished_number, channel.reads.time_spent, 
            channel.reads.number_bps_requested, channel.reads.number_bps_read, channel.reads.number_bps_rejected, channel.reads.fin_number_stop_receiving, 
            channel.reads.fin_number_rejected, channel.reads.number_rejected_missed, channel.reads.number_stop_receiving_missed, 
            channel.no_reads_left.finished_number, channel.no_reads_left.time_spent
        ) for channel in channel_stats], columns=[
            'channel',
            'short_gaps_finished_number', 'short_gaps_time_spent', 'long_gaps_finished_number', 'long_gaps_time_spent', 'unblock_delays_finished_number', 'unblock_delays_time_spent', 
            'mux_scans_finished_number', 'mux_scans_time_spent', 'channel_broken_finished_number', 'channel_broken_time_spent', 
            'reads_finished_number', 'reads_time_spent', 'reads_number_bps_requested', 'reads_number_bps_read', 'reads_number_bps_rejected', 
            'reads_fin_number_stop_receiving', 'reads_fin_number_rejected', 'reads_number_rejected_missed', 'reads_number_stop_receiving_missed', 
            'no_reads_left_finished_number', 'no_reads_left_time_spent'
        ]
    )
    df.sort_values("reads_finished_number", inplace=True)
    
    return df

def plot_read_stats_per_channel(df, save_dir=None):
    """Plot stats about reads per channel, df coming from channel_stats_to_df"""
    df_reads = df[[col for col in df.columns if col.startswith("reads_")]].copy()
    df_reads["reads_number_bps_notrequested"] = df_reads["reads_number_bps_read"] - df_reads["reads_number_bps_requested"] # silently read, but never requested by ReadUntil

    # plot number of basepairs in each channel, splitting by whether they were requested or not, and rejected
    df_reads1 = df_reads[["reads_number_bps_requested", "reads_number_bps_notrequested", "reads_number_bps_rejected"]]
    df_reads1 = df_reads1.rename({"reads_number_bps_requested": "requested by ReadUntil", "reads_number_bps_notrequested": "read but not requested", "reads_number_bps_rejected": "rejected"}, axis=1)
    df_reads1.reset_index(drop=True, inplace=True)
    # fig, ax1 = plt.subplots(figsize=(4 + len(df)/20, 10))
    # df_reads1.plot.bar(stacked=True, ax=ax1)
    # ax1.set_xlabel("Channel (sorted by number of reads)")
    # ax1.set_ylabel("Number of basepairs")
    # ax1.set_title("Number of basepairs in each channel and rejected number of basepairs")
    axes = df_reads1.plot.bar(
        stacked=True, subplots=True, 
        xlabel="Channel (sorted by number of reads)", 
        ylabel="Number of basepairs",
        title="Number of basepairs in each channel and rejected number of basepairs",
        legend=False,
        figsize=(4 + len(df)/20, 10),
    )
    ax1 = axes[0]
    _disable_x_ticks(axes[2])
    fig = ax1.figure
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"nb_basepairs_recrej_per_channel.{FIGURE_EXT}")

    # plot number of reads set to stop_receiving, rejected in each channel
    df_reads2 = df_reads[["reads_fin_number_stop_receiving", "reads_fin_number_rejected"]]
    df_reads2 = df_reads2.rename({"reads_fin_number_stop_receiving": "stop receiving", "reads_fin_number_rejected": "rejected"}, axis=1)
    fig, ax2 = plt.subplots(figsize=(4 + len(df)/20, 10))
    df_reads2.plot.bar(stacked=True, ax=ax2)
    ax2.set_xlabel("Channel (sorted by number of reads)")
    ax2.set_ylabel("Number of successful actions")
    ax2.set_title("Number of successful actions per channel")
    _disable_x_ticks(ax2)
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"nb_successful_stoprec_rejected_per_channel.{FIGURE_EXT}")

    # plot number of reads missed in each channel
    df_reads3 = df_reads[["reads_number_rejected_missed", "reads_number_stop_receiving_missed"]]
    df_reads3 = df_reads3.rename({"reads_number_rejected_missed": "rejected", "reads_number_stop_receiving_missed": "stop receiving"}, axis=1)
    fig, ax3 = plt.subplots(figsize=(4 + len(df)/20, 10))
    df_reads3.plot.bar(stacked=True, ax=ax3)
    ax3.set_xlabel("Channel (sorted by number of reads)")
    ax3.set_ylabel("Number of missed actions")
    ax3.set_title("Number of missed read actions per channel")
    _disable_x_ticks(ax3)
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"nb_missed_stoprec_rejected_per_channel.{FIGURE_EXT}")
    
    return df_reads1, df_reads2, df_reads3, ax1, ax2, ax3

def plot_number_elems_per_channel(df, save_dir=None):
    """Plot barplot of number of elements per element type for each channel"""
    df_number_elems = df[[col for col in df.columns if col.endswith("_finished_number")]].copy()
    df_number_elems.rename(columns={col: col.replace("_finished_number", "") for col in df_number_elems.columns}, inplace=True)
    df_number_elems = df_number_elems[["reads", "short_gaps", "long_gaps", "unblock_delays", "mux_scans", "channel_broken", "no_reads_left"]] # reorder
    df_number_elems

    fig, ax = plt.subplots(figsize=(4 + len(df)/20, 5))
    df_number_elems.plot.bar(stacked=True, ax=ax)
    ax.set_xlabel("Channel (sorted by number of reads)")
    ax.set_ylabel("Number of elements")
    ax.set_title("Number of elements in each channel")
    _disable_x_ticks(ax)
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"nb_elements_per_channel.{FIGURE_EXT}")
    
    return df_number_elems, ax

def plot_time_spent_per_channel(df, save_dir=None):
    """Plot time spent in each element type per channel"""
    
    df_time_spent = df[[col for col in df.columns if col.endswith("_time_spent")]].copy()
    df_time_spent.rename(columns={col: col.replace("_time_spent", "") for col in df_time_spent.columns}, inplace=True)
    df_time_spent = df_time_spent[["reads", "short_gaps", "long_gaps", "unblock_delays", "mux_scans", "channel_broken", "no_reads_left"]] # reorder
    
    fig, ax = plt.subplots(figsize=(4 + len(df)/20, 5))
    df_time_spent.plot.bar(stacked=True, ax=ax)
    ax.set_xlabel("Channel (sorted by number of reads)")
    ax.set_ylabel("Time spent (s)")
    ax.set_title("Time spent in each element type per channel")
    _disable_x_ticks(ax)
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"time_spent_per_element_per_channel.{FIGURE_EXT}")

    return df_time_spent, ax

def plot_channel_stats(channel_stats: List[ChannelStats], save_dir=None, close_figures=None):
    """
    Per channel, plot read stats, number of elements, time spent per element
    
    Args:
        channel_stats: list of ChannelStats
        save_dir: directory to save figures to
        close_figures: whether to close each figure after saving it (reduces memory usage as this function plots a lot)
        
    Returns:
        data frame with one line per channel
    """
    if close_figures is None:
        close_figures = save_dir is not None
        
    def close_fig(fig):
        if close_figures:
            plt.close(fig)
            
    df = channel_stats_to_df(channel_stats)
    ax1, ax2, ax3 = plot_read_stats_per_channel(df, save_dir=save_dir)[3:]; logger.debug("Created 3 plots"); close_fig(ax1.figure); close_fig(ax2.figure); close_fig(ax3.figure)
    ax = plot_number_elems_per_channel(df, save_dir=save_dir)[1]; logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_time_spent_per_channel(df, save_dir=save_dir)[1]; logger.debug("Created 1 plot"); close_fig(ax.figure)
    
    return df