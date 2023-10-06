"""
Analyze and plot a sequencing summary file from the simulator

If save_dir is provided, the plots are closed to save memory.
"""

import argparse
from collections import Counter
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import warnings
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import tqdm
from intervaltree import Interval, IntervalTree

from simreaduntil.seqsum_tools.seqsum_preprocessing import add_previous_gap_duration, compute_long_gap_threshold, get_gaps_single_channel, sort_and_clean_seqsum_df
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, setup_logger_simple
from simreaduntil.shared_utils.merge_axes import save_fig_and_pickle
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId
from simreaduntil.seqsum_tools.coverage_tracker import NanoSimCoverageTracker, PafCoverageTracker
from simreaduntil.shared_utils.plotting import make_tight_layout, rotate_xticks
from simreaduntil.shared_utils.plotting import FIGURE_EXT
from simreaduntil.shared_utils.utils import is_sorted, print_args
from simreaduntil.shared_utils.debugging_helpers import is_test_mode
from simreaduntil.shared_utils.plotting import _disable_x_ticks

GROUP_COLUMN = "chrom"
"""column to group by, e.g. chromosome"""

logger = setup_logger_simple(__name__)
"""module logger"""

def set_plots_groupby_column(group_column):
    """
    Set the column to group by, for this module
    
    This method has side effects.
    
    Args:
        group_column: column to group by, e.g. chromosome
    """
    global GROUP_COLUMN
    GROUP_COLUMN = group_column
    
def add_group_and_reflen_from_nanosim_id(seqsum_df, group_column=GROUP_COLUMN):
    """
    Parse the chromosome (as group) and reference length from the NanoSim read ids, in-place
    
    Args:
        seqsum_df: seqsummary dataframe
        group_column: column to group by, e.g. chromosome
    """
    group_column = group_column or GROUP_COLUMN
    
    # assigning two columns is very slow: takes 14s vs 0.4s below
    # def get_chrom_and_ref_len(read_id):
    #     nanosim_id = NanoSimId.from_str(read_id)
    #     return pd.Series([nanosim_id.chrom, nanosim_id.ref_len])
    # seqsum_df[[group_column, "nb_ref_bps"]] = seqsum_df["read_id"].apply(get_chrom_and_ref_len)
    
    def get_chrom(read_id):
        nanosim_id = NanoSimId.from_str(read_id)
        return nanosim_id.chrom

    def get_ref_len(read_id):
        nanosim_id = NanoSimId.from_str(read_id)
        return nanosim_id.ref_len

    seqsum_df[group_column] = seqsum_df["read_id"].apply(get_chrom)
    seqsum_df["nb_ref_bps"] = seqsum_df["read_id"].apply(get_ref_len)
    
def add_group_and_reflen_from_paf(seqsum_df, paf_file, group_column=None):
    """
    Adds the group (to group by) as well as the reference length (nb_ref_bps) from the PAF file
    
    If the read ids are NanoSim ids, you should not call this function if the sequencing summary file comes from the simulator since it should already contain it.
    ReadFish's sequencing summary file usually contains a subset of the ids because it may not see all reads, therefore use minimap2 to map the reads to the reference.
    Reads that are in the sequencing summary file but not in the PAF file receive group_column 'unmapped' and 'nb_ref_bps' is set to 'sequence_length_template'.
    
    Args:
        seqsum_df: seqsummary dataframe
        paf_file: PAF file without duplicated ids (using --secondary=no in minimap2)
        group_column: column to group by, e.g. chromosome
        
    Returns:
        seqsum_df with additional columns 'group_column' and 'nb_ref_bps'; must not initially be present
    """
    group_column = group_column or GROUP_COLUMN
            
    assert (group_column not in seqsum_df.columns) and ("nb_ref_bps" not in seqsum_df.columns), f"columns {group_column} and nb_ref_bps must not be present"

    # print warning if the first read_id is a valid NanoSim id
    if NanoSimId.is_valid(seqsum_df["read_id"].iloc[0]):
        warnings.warn("The read ids are NanoSim ids, use 'add_group_and_reflen_from_nanosim_id' instead. Note that the NanoSim ids of rejected reads are modified by the simulator.")

    # column 6 (1-based): target sequence name, 8: target start, 9: target end (it seems excluding this index)
    logger.debug(f"Reading PAF file '{paf_file}'")
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=[0, 5, 7, 8], names=["read_id", group_column, "target_start", "target_end"])
    logger.debug("Finished reading PAF file")
    
    assert paf_df["read_id"].is_unique, "read ids in PAF file must be unique, probably took PAF file from ReadFish which logs each time a chunk maps" # see https://github.com/lh3/minimap2/blob/master/FAQ.md to filter out chimeric alignments
    paf_df["nb_ref_bps"] = paf_df["target_end"] - paf_df["target_start"] # no +1 offset, interval [a, b)
    paf_df.drop(["target_start", "target_end"], axis=1, inplace=True)

    # left join based on read_id
    seqsum_df = seqsum_df.set_index("read_id", append=True)
    paf_df.set_index("read_id", inplace=True)
    seqsum_df = seqsum_df.join(paf_df, how="left", on="read_id")
    fraction_unmapped = seqsum_df[group_column].isna().sum() / len(seqsum_df)
    logger.info(f"Did not map {fraction_unmapped:.2%} of reads")
    na_rows = seqsum_df[group_column].isna() # reads not mapped to a genome
    seqsum_df[group_column].fillna("unmapped", inplace=True)
    seqsum_df.loc[na_rows, "nb_ref_bps"] = seqsum_df.loc[na_rows, "sequence_length_template"] # approximate ref length (ignoring indels)
    seqsum_df["nb_ref_bps"] = seqsum_df["nb_ref_bps"].astype(int)
    return seqsum_df.reset_index(level="read_id")
        
def seqsum_add_cols_for_plotting_selseq_performance(seqsum_df, group_column=None) -> pd.DataFrame:
    """
    Add columns for plotting, e.g. cumulative number of rejections, per group
    
    The data frame is sorted by end time to compute the cumulative number of rejections until that tim, per group.
    
    Args:
        seqsum_df: seqsummary dataframe
        group_column: column to group by, e.g. chromosome
    
    Returns:
        seqsum_df with additional columns to plot
    """
    group_column = group_column or GROUP_COLUMN
    
    if "end_reason" not in seqsum_df.columns:
        seqsum_df["end_reason"] = "unknown"
    seqsum_df["end_reason"] = seqsum_df["end_reason"].astype("category")
        
    seqsum_df["is_user_rejection"] = seqsum_df["end_reason"] == "data_service_unblock_mux_change" # user rejections only, rejections due to mux scan may also happen
    seqsum_df["is_full_read"] = seqsum_df["end_reason"] == "signal_positive" # number of full reads
    seqsum_df["end_time"] = seqsum_df["start_time"] + seqsum_df["duration"]
    
    assert "nb_ref_bps" in seqsum_df.columns, f"column nb_ref_bps must be present"
    assert group_column in seqsum_df.columns, f"column '{group_column}' must be present"
    seqsum_df[group_column] = seqsum_df[group_column].astype("category") # seaborn hue arg requires category (otherwise sometimes does not show all categories, even if non-empty)
        
    # required for cumulative numbers
    seqsum_df = seqsum_df.sort_values("end_time")
        
    if "nb_ref_bps_full" in seqsum_df.columns:
        seqsum_df["nb_rejectedbps"] = seqsum_df["nb_ref_bps_full"] - seqsum_df["nb_ref_bps"]

    # fields starting with cum_nb_ are cumulative
    seqsum_df[f"cum_nb_reads_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True).cumcount() + 1
    seqsum_df[f"cum_nb_rejections_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["is_user_rejection"].cumsum()
    seqsum_df[f"cum_nb_full_reads_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["is_full_read"].cumsum()
    if "nb_ref_bps_full" in seqsum_df.columns:
        seqsum_df[f"cum_nb_rejectedbps_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["nb_rejectedbps"].cumsum()
    seqsum_df[f"cum_nb_ref_bps_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["nb_ref_bps"].cumsum()
    seqsum_df[f"cum_nb_seq_bps_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["sequence_length_template"].cumsum()
    if "stopped_receiving" in seqsum_df.columns:
        seqsum_df[f"cum_nb_stopped_receiving_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["stopped_receiving"].cumsum()
    if "never_requested" in seqsum_df.columns:
        seqsum_df[f"cum_nb_never_requested_per_{group_column}"] = seqsum_df.groupby(group_column, observed=True)["never_requested"].cumsum()

    if "nb_ref_bps_full" in seqsum_df.columns:
        assert all(~seqsum_df["is_user_rejection"] | seqsum_df["nb_rejectedbps"] > 0) # user rejected => nb_rejectedbps > 0
    
    return seqsum_df

# todo2: test group_to_chroms arg
def compute_coverage_per_group_df(seqsum_df, cov_tracker: NanoSimCoverageTracker, cov_thresholds=[1, 2, 3, 4, 5, 6], coverage_every=1, chrom_column="chrom", group_to_chroms=None):
    """
    Compute coverage per group at regular intervals
    
    This function assumes that group_column corresponds to the chromosome.
    It adds a column named "group".
    
    Args:
        seqsum_df: seqsummary dataframe
        cov_tracker: initial cov tracker to which reads are added
        cov_thresholds: target coverages to compute
        coverage_every: compute coverage every coverage_every reads
        chrom_column: column containing chromosome info
        group_to_chroms: chromosomes per group; if None, separate group per chromosome except "unmapped"
        
    Returns:
        dataframe with coverage per group at regular intervals, with "group" column
    """
    assert chrom_column in seqsum_df.columns
    assert is_sorted(seqsum_df["end_time"].values)
    
    if group_to_chroms is None:
        group_to_chroms = {x: [x] for x in seqsum_df[chrom_column].unique() if x != "unmapped"}
    logger.debug(f"Computing coverage per group {group_to_chroms} every {coverage_every} reads")
    cov_data = []
    for (i, (read_id, end_time)) in enumerate(tqdm.tqdm(zip(seqsum_df["read_id"], seqsum_df["end_time"]), desc="Adding reads to coverage tracker", total=len(seqsum_df))):
        cov_tracker.add_read(read_id)
        if i % coverage_every == 0:
            for (group_name, chroms) in group_to_chroms.items():
                cov_data.append([end_time, group_name] + list(cov_tracker.get_fraction_cov_atleast(cov, chroms=chroms) for cov in cov_thresholds))
    
    logger.debug("Done computing coverage")
    return pd.DataFrame(cov_data, columns=["end_time", "group"] + [f"fraction_cov{cov}" for cov in cov_thresholds])

def plot_number_reads_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of reads"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_reads_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of reads")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_reads_per_{group_column}.{FIGURE_EXT}")
        
    return ax


def plot_number_full_reads_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of full reads"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_full_reads_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of full reads")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_full_reads_per_{group_column}.{FIGURE_EXT}")
        
    return ax

def plot_number_rejections_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of user rejections"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_rejections_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of user rejections")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_rejections_per_{group_column}.{FIGURE_EXT}")
        
    return ax

def plot_number_stop_receiving_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of reads successfully set to stop receiving"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_stopped_receiving_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of reads set to stop receiving")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_stopped_receiving_per_{group_column}.{FIGURE_EXT}")
        
    return ax

def plot_number_never_requested_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of reads that were never requested by readuntil"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_never_requested_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of reads that were never requested by ReadUntil")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_never_requested_per_{group_column}.{FIGURE_EXT}")
        
    return ax

def plot_fraction_never_requested_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot fraction of reads that were never requested by readuntil"""
    group_column = group_column or GROUP_COLUMN
    
    seqsum_df[f"cum_frac_never_requested_per_{group_column}"] = seqsum_df[f"cum_nb_never_requested_per_{group_column}"] / seqsum_df[f"cum_nb_reads_per_{group_column}"]
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_frac_never_requested_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Fraction of reads that were never requested by ReadUntil")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_never_requested_per_{group_column}.{FIGURE_EXT}")
        
    return ax

def plot_number_rejectedbps_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of basepairs rejected (with respect to full read)"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_rejectedbps_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of rejected basepairs") # not only by ReadUntil reject, but also due to mux scan
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_rejectedbps_per_{group_column}.{FIGURE_EXT}")
        
    return ax
    
def plot_number_sequenced_bps_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of basepairs over time"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_seq_bps_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of sequenced basepairs")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_seq_bps_per_{group_column}.{FIGURE_EXT}")
    
    return ax
    
def plot_number_ref_bps_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of basepairs (reference bps) over time"""
    group_column = group_column or GROUP_COLUMN
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(seqsum_df, x="end_time", y=f"cum_nb_ref_bps_per_{group_column}", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Number of sequenced basepairs (reference basepairs)") # number of bps on reference, not those read (differs due to indels)
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"cum_nb_ref_bps_per_{group_column}.{FIGURE_EXT}")
        
    return ax
    
def plot_read_duration_per_group(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot read duration over time"""
    group_column = group_column or GROUP_COLUMN
    
    # fig, ax = plt.subplots(figsize=(10, 5))
    # sns.scatterplot(seqsum_df, x="end_time", y="duration", hue=group_column, ax=ax, s=1)
    
    sns.lmplot(seqsum_df, x="end_time", y="duration", hue=group_column, scatter_kws={"s": 1})
    ax = plt.gca()
    fig = plt.gcf()
    
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel("Read duration (s)")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"read_duration.{FIGURE_EXT}")
        
    return ax
        
def plot_read_length_per_group_per_decisiontype(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot number of basepairs in read over time"""
    group_column = group_column or GROUP_COLUMN
    
    valid_subset = {"signal_positive", "data_service_unblock_mux_change"}
    if not set(seqsum_df["end_reason"].unique()).issubset(valid_subset):
        # real seqsum file has values: signal_positive, data_service_unblock_mux_change, signal_negative, unblock_mux_change, mux_change
        warnings.warn(f"""Only plotting full reads (signal_positive) and user rejected reads (data_service_unblock_mux_change). There are reads with other end reasons as well: {set(seqsum_df["end_reason"].unique()).difference(valid_subset)}.""")
    
    for (extra_text, sub_df) in [("rejected reads", seqsum_df[seqsum_df["is_user_rejection"]]), ("full reads", seqsum_df[seqsum_df["is_full_read"]])]:
        
        # fig, ax = plt.subplots(figsize=(10, 5))    
        # sns.scatterplot(sub_df, x="end_time", y="sequence_length_template", hue=group_column, ax=ax, s=1)
        
        sns.lmplot(sub_df, x="end_time", y="sequence_length_template", hue=group_column, scatter_kws={"s": 1})
        ax = plt.gca()
        fig = plt.gcf()
        
        ax.set_xlabel("Read end time (s)")
        ax.set_ylabel(f"Read length (bp) of {extra_text}")
        make_tight_layout(fig)
        if save_dir is not None:
            save_fig_and_pickle(fig, save_dir / f"read_length_{extra_text.split(' ')[0]}.{FIGURE_EXT}")
            
    return ax
    
def plot_fraction_read_per_group_per_rejected(seqsum_df, save_dir: Optional[Path]=None, group_column=None):
    """Plot fraction of basepairs read over time, for rejected reads"""
    group_column = group_column or GROUP_COLUMN
    
    seqsum_df["fraction_read"] = seqsum_df["nb_ref_bps"] / seqsum_df["nb_ref_bps_full"]
    fig, ax = plt.subplots(figsize=(10, 5))
    df1 = seqsum_df[seqsum_df["is_user_rejection"]]
    if len(df1) > 0:
        # seaborn bug when data frame has 0 rows and specifying hue does not work with seaborn0.13.0, but with seaborn0.12.2
        sns.lineplot(df1, x="end_time", y="fraction_read", hue=group_column, ax=ax)
    ax.set_xlabel("Read end time (s)")
    ax.set_ylabel(f"Fraction of read until rejection, for rejected reads")
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"fraction_until_rejection.{FIGURE_EXT}")
        
    return ax
    
def plot_number_channels_per_group_over_time(seqsum_df, timepoints=None, group_column=None, save_dir=None):
    """
    Plot how many channels are reading a read from a group over time, creates one normalized and one unnormalized plot
    """
    group_column = group_column or GROUP_COLUMN

    # add something to end since upper bound is not inclusive
    reads_intervaltree = IntervalTree(Interval(start_time, end_time, group) for (start_time, end_time, group) in zip(seqsum_df["start_time"], seqsum_df["end_time"] + 1e-8, seqsum_df[group_column]))
    if timepoints is None:
        timepoints = np.linspace(seqsum_df["start_time"].min(), seqsum_df["end_time"].max(), 500)

    def count_per_group(read_objs):
        return Counter(o.data for o in read_objs)
        # res = Counter(o.data for o in read_objs)
        # res["other"] = n_channels - sum(res.values())
        # return res

    frac_per_group_over_time = np.array([count_per_group(reads_intervaltree[t]) for t in timepoints])
    
    def create_plot(ax, normalize):
        # normalize: whether to normalize by number of currently reading channels
        df_proportions = pd.DataFrame.from_records(frac_per_group_over_time) # wide data-frame with one row per timepoint
        if normalize:
            # if only NaNs in row, sum is 0 for that row
            # divide row-wise
            df_proportions = df_proportions.divide(df_proportions.sum(axis=1, skipna=True), axis=0)
        df_proportions.fillna(0, inplace=True)
        df_proportions["time"] = timepoints

        df_proportions.set_index("time", inplace=True)
        df_proportions.plot.area(ax=ax)
        # multiple argument not available for sns.barplot
        # df1 = pd.melt(df_proportions, id_vars=["time"], var_name=group_column, value_name="nb_reads")
        # sns.barplot(df1, x="time", y="nb_reads", hue=group_column, ax=ax, multiple="stack")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(f"{'Proportion' if normalize else 'Number'} of active channels\nreading from {group_column}")
        ax.set_title(f"""{'Proportion' if normalize else 'Number'} of active channels reading from {group_column} over time (stacked), {seqsum_df["channel"].nunique()} active channels""")
        # df_proportions.plot.bar(stacked=True, rot=90, ax=ax)
        
        return ax
    
    fig, (ax1, ax2) = plt.subplots(figsize=(12, 8), nrows=2, sharex=True)
    create_plot(ax1, normalize=False)
    create_plot(ax2, normalize=True)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"proportion_channels_per_{group_column}_over_time.{FIGURE_EXT}")
        
    return fig

def plot_channel_occupation_fraction_over_time(seqsum_df, timepoints=None, mux_scan_interval=None, save_dir=None):
    """
    Plot channel occupation (active percentage) over time
    
    A channel is active until it has produced its last read.
        
    Args:
        seqsum_df: dataframe with columns "start_time", "end_time", "channel"
        timepoints: timepoints to plot at; if None, a linspace between the min and max start_time is used
        mux_scan_interval: interval at which the mux scan is performed (in seconds); if None, no vertical lines are plotted
        save_dir: directory to save plot to; if None, not saving plot
        
    Returns:
        axis object
    """
    
    n_channels = seqsum_df["channel"].nunique()
    
    
    # we don't need to store the channel as data in the IntervalTree, so we use .from_tuples()
    read_interval_tree = IntervalTree.from_tuples(zip(seqsum_df["start_time"], seqsum_df["end_time"] + 1e-8)) # add something to end since upper bound is not inclusive

    # plot fraction of active channels over time
    if timepoints is None:
        timepoints = np.linspace(seqsum_df["start_time"].min(), seqsum_df["end_time"].max(), 500)
    
    n_channels_active_per_time = np.array([len(read_interval_tree[t]) for t in timepoints])
    
    # find total number of basepairs read at the timepoints
    # cum_nb_bps = seqsum_df["sequence_length_template"].cumsum()
    # cum_nb_bps_at_timepoints = cum_nb_bps.iloc[np.maximum(0, np.searchsorted(seqsum_df["end_time"], timepoints, side="right") - 1)]

    fig, ax = plt.subplots(figsize=(20, 10))
    # ax.plot(timepoints, n_channels_active_per_time / n_channels * 100) # note: this is the fraction among the active channels, not including broken channels (i.e. those never producing any read)
    ax.plot(timepoints, n_channels_active_per_time)
    # ax.plot(cum_nb_bps_at_timepoints, n_channels_active_per_time / n_channels * 100)
    
    if mux_scan_interval is not None:
        ax.vlines(x=[mux_scan_interval*i for i in range(int(np.ceil(max(timepoints)/mux_scan_interval)))], ymin=0, ymax=1, label="mux scans every 90 mins", color="red", linestyle="--", transform=ax.get_xaxis_transform())
    
    # ax.set_xlim(0)
    ax.set_ylim(0)
    ax.set_xlabel("Time (s)")
    # ax.set_xlabel("Number of basepairs across all channels")
    # ax.set_ylabel("Fraction of active channels (%)")
    # ax.set_title(f"Fraction of active channels over time ({n_channels} active channels)")
    ax.set_ylabel("Number of reading channels")
    ax.set_title(f"Number of reading channels over time ({n_channels} active channels)")
    ax.legend()
    make_tight_layout(fig)
    
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"channel_occupation_fraction_over_time.{FIGURE_EXT}")
    
    return ax

def keep_largest_gaps_only(segments, k):
    """
    Keep the largest k gaps (or less if there are less gaps) between segments, remove the smaller gaps by merging the corresponding segments.
    
    Args:
        segments: list of segments [[(x_start, y), (x_end, y)], ..]
            assumes y is the same for all entries
        k: maximum number of gaps to keep
    """
    assert k >= 0
    if len(segments) <= k+1:
        return segments
    # at least two segments
    
    y = segments[0][0][1]
    segments_x = np.array([[seg[0][0], seg[1][0]] for seg in segments])
    
    gap_starts = segments_x[:-1, 1]
    gap_ends = segments_x[1:, 0]
    gap_indices_to_keep = np.argsort(gap_ends - gap_starts)[len(gap_ends)-k:] # just '-k' not working when k=0
    gap_indices_to_keep.sort()
    new_segments_x = np.stack((
        np.concatenate(([segments_x[0][0]], gap_ends[gap_indices_to_keep])),
        np.concatenate((gap_starts[gap_indices_to_keep], [segments_x[-1][1]]))
    ), axis=1)
    
    return np.array([[(seg[0], y), (seg[1], y)] for seg in new_segments_x])

def plot_channels_over_time(seqsum_df, long_gap_threshold=None, figsize=None, save_dir=None, max_num_gaps_per_channel=50) -> plt.Axes:
    """
    Plot reads and gaps over time for all channels
    
    Note: It is better to use the argument "max_num_gaps_per_channel" than setting it to None and storing the image as a 
    low-resolution png (to save space) since the plotting still is done before which is time-consuming.
    
    Args:
        seqsum_df: dataframe with columns "start_time", "end_time", "channel"
        long_gap_threshold: threshold for long gaps (in seconds); if None, computed from data (requires "prev_gap_duration" column)
        figsize: figure size
        save_dir: directory to save plot to; if None, not saving plot
        max_num_gaps_per_channel: maximum number of gaps to plot per channel (to reduce number of points to plot); if None, all gaps are plotted
    """
    
    n_channels = seqsum_df["channel"].nunique()
    
    if figsize is None:
        figsize = 40, 1 + 7/60 * n_channels
    fig, ax = plt.subplots(figsize=figsize)
    # fig, ax = plt.subplots(figsize=(200, 30))
    
    if long_gap_threshold is None:
        long_gap_threshold = compute_long_gap_threshold(seqsum_df["prev_gap_duration"])
    
    # plot reads in red at position y_pos (typically df corresponds to a single channel)
    def plot_read_segments(df_single, y_pos):
        segments =[((start, y_pos), (end, y_pos)) for (start, end) in zip(df_single["start_time"], df_single["end_time"])]
        # lines.append([(elem.t_start, y_pos + offset), (elem.t_end, y_pos + offset)])
        if max_num_gaps_per_channel is not None:
            segments = keep_largest_gaps_only(segments, k=max_num_gaps_per_channel)
        line_collection = LineCollection(segments, linewidths=5, colors="red")
        ax.add_collection(line_collection)
        
    # plot gaps coming after reads, in blue, linewidth 2
    def plot_gap_segments(gap_starts, gap_durations, y_pos):
        if max_num_gaps_per_channel is not None:
            gap_indices_to_keep = np.argsort(gap_durations)[-max_num_gaps_per_channel:]
            gap_starts, gap_durations = gap_starts[gap_indices_to_keep], gap_durations[gap_indices_to_keep]
        
        segments = [((start, y_pos), (end, y_pos)) for (start, end) in zip(gap_starts, gap_starts+gap_durations)]
        # lines.append([(elem.t_start, y_pos + offset), (elem.t_end, y_pos + offset)])
        line_collection = LineCollection(segments, linewidths=2, colors="blue")
        ax.add_collection(line_collection)
        
    df_channel_groups = seqsum_df.groupby("channel", observed=True)
    for (channel, df_single) in tqdm.tqdm(df_channel_groups, total=len(df_channel_groups), desc="Plotting channels over time"):
        y_pos = channel
        
        # plot reads that pass filtering separately from those that don't
        plot_read_segments(df_single[df_single["passes_filtering"]], y_pos-0.05)
        plot_read_segments(df_single[~df_single["passes_filtering"]], y_pos+0.05)
        channel_gap_durations, channel_gap_starts = get_gaps_single_channel(df_single, seq_start_time=0, check_no_overlap=False)
        plot_gap_segments(channel_gap_starts[channel_gap_durations > long_gap_threshold], channel_gap_durations[channel_gap_durations > long_gap_threshold], y_pos=y_pos+0.02)
        
        # plot a line from first to last time point
        first_time = df_single["start_time"].min()
        last_time = df_single["end_time"].max()
        # ax.plot([first_time, last_time], [y_pos+0.05, y_pos+0.05], color="gray", alpha=0.5, linewidth=5)
        ax.fill_between([first_time, last_time], [y_pos-0.05, y_pos-0.05], [y_pos+0.05, y_pos+0.05], color="gray", alpha=0.5)
        # plot a dot at first and last time point
        ax.plot(first_time, y_pos, marker=".", markersize=10, color="green", alpha=0.5)
        ax.plot(last_time, y_pos, marker=".", markersize=10, color="black", alpha=0.5)
        
        # if channel > 2:
        #     break

    # create legend for reads and gaps
    existing_point = (first_time, y_pos)
    legend_elements = [
        Line2D(existing_point, existing_point, color='red', lw=2, label='read'),
        Line2D(existing_point, existing_point, color='blue', lw=2, label='long gap'),
        Line2D(existing_point, existing_point, marker=".", markersize=10, color="green", linestyle="None", label='start time'),
        Line2D(existing_point, existing_point, marker=".", markersize=10, color="black", linestyle="None", label='end time'),
    ]
    ax.legend(handles=legend_elements, loc='center right')
    
    ax.autoscale() # Manually adding artists doesn't rescale the plot, so we need to autoscale
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Channel")
    ax.set_title(f"Reads over time ({n_channels} active channels)")
    # make yticks at integer values only
    # ax.set_yticks(np.arange(1, n_channels+1))
    make_tight_layout(fig)

    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"channels_over_time.{FIGURE_EXT}")
        
    return ax

def plot_mux_over_time(seqsum_df, figsize=None, save_dir=None, max_points_per_channel=20):
    """
    Plots the mux over time per channel
    
    Each channel is plotted at a y position and the offset indicates the mux.
    
    Args:
        seqsum_df: dataframe with columns "start_time", "end_time", "channel", "mux"
        figsize: figure size
        save_dir: directory to save plot to; if None, not saving plot
        max_points_per_channel: maximum number of points to plot per channel by sampling (to reduce plotting time); if None, do not reduce
    
    Returns:
        axis object
    """
    n_channels = seqsum_df["channel"].nunique()
    
    if figsize is None:
        figsize = 40, 1 + 7/60 * n_channels

    fig, ax = plt.subplots(figsize=figsize)

    df_channel_groups = seqsum_df.groupby("channel", observed=True)
    for (channel, df_single) in tqdm.tqdm(df_channel_groups, total=len(df_channel_groups), desc="Processing channels for mux plot"):
        y_pos = channel
        
        def get_ypos_for_mux(mux):
            return y_pos + (mux - 2.5) * 0.3
        
        if max_points_per_channel is not None:
            df_single = df_single.sample(min(len(df_single), max_points_per_channel))

        # mux is from 1 to 4
        ax.plot(df_single["start_time"], get_ypos_for_mux(df_single["mux"]), color="blue", marker=".", markersize=3)
        
        # # plot grey lines for different muxes
        # for mux in range(1, 5):
        #     ax.axhline(y=get_ypos_for_mux(mux), color="gray", alpha=0.25 + 0.75*(mux/4))
        
        # if channel >= 2:
        #     break
        
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Channel + mux")
    ax.set_title("Mux (well in channel) over time for the channel")
    make_tight_layout(fig)
    
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"mux_over_time.{FIGURE_EXT}")
        
    return ax

def plot_read_stats_by_channel_hists(seqsum_df, nbins=30, hue_col=None, save_dir=None):
    """
    Plot statistics about the reads, e.g. number of reads per channel, number of reads per second etc.
    """
    # assume reads are non-overlapping, so can just take the sum
    df_read_stats = seqsum_df.groupby("channel", observed=True).apply(lambda df: pd.Series({"reading_time": df["duration"].sum(), "time_sequencing": df["end_time"].max() - df["start_time"].min(), "num_reads": df.shape[0], "num_bps": df["sequence_length_template"].sum()}))
    df_read_stats

    df_read_stats["reads_per_second"] = df_read_stats["num_reads"] / df_read_stats["time_sequencing"] # only counting the active time
    df_read_stats["perc_time_reading"] = df_read_stats["reading_time"] / df_read_stats["time_sequencing"] * 100
    def add_mean_median(vals, ax):
        ax.axvline(x=np.mean(vals), linestyle="dotted", color="red", label="mean")
        ax.axvline(x=np.median(vals), linestyle="dotted", color="black", label="median")
    
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5, figsize=(22, 5))
    
    sns.histplot(df_read_stats, x="time_sequencing", hue=hue_col, bins=nbins, ax=ax1)
    add_mean_median(df_read_stats["time_sequencing"], ax=ax1)
    ax1.set_xlabel("Time spent sequencing (s)")
    ax1.set_ylabel("Number of channels")
    ax1.legend()
    rotate_xticks(ax1, rotation=45)
    
    sns.histplot(df_read_stats, x="num_bps", hue=hue_col, bins=nbins, ax=ax2)
    add_mean_median(df_read_stats["num_bps"], ax=ax2)
    ax2.set_xlabel("Number of sequenced basepairs")
    ax2.set_ylabel("Number of channels")
    ax2.legend()
    rotate_xticks(ax2, rotation=45)
    
    sns.histplot(df_read_stats, x="num_reads", hue=hue_col, bins=nbins, ax=ax3)
    ax3.set_xlabel("Number of reads in channel")
    ax3.set_ylabel("Number of channels")
    add_mean_median(df_read_stats["num_reads"], ax=ax3)
    ax3.legend()
    
    sns.histplot(df_read_stats, x="perc_time_reading", hue=hue_col, bins=nbins, ax=ax4)
    add_mean_median(df_read_stats["perc_time_reading"], ax=ax4)
    ax4.set_xlabel("Percentage of time spent reading")
    ax4.set_ylabel("Number of channels")
    ax4.legend()
    
    sns.histplot(df_read_stats, x="reads_per_second", hue=hue_col, bins=nbins, ax=ax5)
    ax5.set_xlabel("Reads per second while active")
    ax5.set_ylabel("Number of channels")
    add_mean_median(df_read_stats["reads_per_second"], ax=ax5)
    ax5.legend()
        
    fig.suptitle("Read statistics with histograms over channels")
    make_tight_layout(fig)
    
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"read_stats_by_channel.{FIGURE_EXT}")
    
    return fig, df_read_stats

def plot_fraction_states_per_channel(seqsum_df, long_gap_threshold=None, end_time=None, save_dir=None):
    """
    Plot time spent reading, in long gaps and short gaps per channel.
    
    Mux scans should be removed beforehand.
    
    Args:
        seqsum_df: sequencing summary dataframe
        long_gap_threshold: threshold for long gaps; if None, it is computed from the data
        end_time: end time of the sequencing run; if None, it is computed from the data
        save_dir: if not None, save the figure to this directory
    """
    if long_gap_threshold is None:
        long_gap_threshold = compute_long_gap_threshold(seqsum_df["prev_gap_duration"])
        
    if end_time is None:
        end_time = seqsum_df["end_time"].max()
    else:
        assert end_time >= seqsum_df["end_time"].max()
    
    def get_time_fractions(seqsum_df):
        gaps = np.concatenate((seqsum_df["start_time"], [end_time])) - np.concatenate(([0], seqsum_df["end_time"].values))
        return pd.Series({
            "reading_time": seqsum_df["duration"].sum(),
            "short_gap_duration": gaps[gaps <= long_gap_threshold].sum(),
            "long_gap_duration": gaps[gaps > long_gap_threshold].sum(),
        })
    
    fractions_per_channel_df = seqsum_df.groupby("channel", observed=True).apply(get_time_fractions)
    fractions_per_channel_df.sort_values("reading_time", ascending=False, inplace=True)
    ax = fractions_per_channel_df.plot.bar(
        stacked=True,
        xlabel="Channel",
        ylabel="Time spent (s)",
        title=f"Time spent per channel (long_gap_threshold: {long_gap_threshold:.2e}s)",
        figsize=(4 + len(fractions_per_channel_df)/20, 5),
    )
    _disable_x_ticks(ax)
    make_tight_layout(ax.figure)
    
    if save_dir is not None:
        save_fig_and_pickle(ax.figure, save_dir / f"time_spent_per_channel.{FIGURE_EXT}")
    
    return ax

def plot_read_end_reason_hist(seqsum_df, save_dir=None):
    """
    Plot histogram over read end reasons
    """
    fig, ax = plt.subplots()
    sns.histplot(seqsum_df["end_reason"], ax=ax)
    rotate_xticks(ax, rotation=45)
    make_tight_layout(fig)
    
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"end_reason_hist.{FIGURE_EXT}")
    
    return ax

def plot_processed_seqsum(seqsum_df, save_dir: Optional[Path]=None, group_column=None, close_figures: Optional[bool]=None):
    """
    Plot a bunch of stuff from the processed seqsum_df, subsampling when sensible
    
    Use ignore_tight_layout_warning to silence the warnings from matplotlib about the tight_layout.
    
    Args:
        seqsum_df: processed seqsum_df
        save_dir: directory to save plots to; if None, not saving them
        group_column: column to group by, e.g. chromosome
        close_figures: whether to close each figure after saving it (reduces memory usage as this function plots a lot)
    """
    group_column = group_column or GROUP_COLUMN
    if close_figures is None:
        close_figures = save_dir is not None
        
    def close_fig(fig):
        if close_figures:
            plt.close(fig)
           
    # # compute instantaneous rates
    # # take difference (x[i+step] - x[i-step]) while keeping array size the same by extending the array on both sides by the step
    # take_diff = lambda x, step=3: np.concatenate((x[step:] - x[:-step], [np.NaN]*step))
    # seqsum_df[f"instrate_reads_per_{group_column}"] = take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100
    # seqsum_df[f"instrate_rejections_per_{group_column}"] = take_diff(seqsum_df[f"cum_nb_rejections_per_{group_column}"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100
    # seqsum_df[f"instrate_rejectedbps_per_{group_column}"] = take_diff(seqsum_df[f"cum_nb_rejectedbps_per_{group_column}"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100
    # seqsum_df[f"instrate_stopped_receiving_per_{group_column}"] = take_diff(seqsum_df[f"cum_nb_stopped_receiving_per_{group_column}"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100
    # seqsum_df[f"instrate_never_requested_per_{group_column}"] = take_diff(seqsum_df[f"cum_nb_never_requested_per_{group_column}"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100
    # seqsum_df[f"instrate_ref_bps"] = take_diff(seqsum_df["cum_nb_ref_bps"].values) / take_diff(seqsum_df[f"cum_nb_reads_per_{group_column}"].values) * 100

    # subsampling to reduce plotting time, but per group to avoid late first read for a particular group that does not occur very often (e.g. chr20,21)
    def sample_group(group):
        return group.sample(n=min(500, len(group)))
    subsampled_seqsum_df = seqsum_df.groupby(group_column, group_keys=False, observed=True).apply(sample_group)
    # subsampled_seqsum_df = seqsum_df.sample(n=min(len(seqsum_df), 500))
    # subsampled_seqsum_df = seqsum_df.reset_index().sample(n=500).sort_index() # keep order, not required here
    
    ax = plot_number_reads_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_number_full_reads_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_number_rejections_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    
    if f"cum_nb_stopped_receiving_per_{group_column}" in seqsum_df.columns:
        ax = plot_number_stop_receiving_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    if f"cum_nb_never_requested_per_{group_column}" in seqsum_df.columns:
        ax = plot_number_never_requested_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
        ax = plot_fraction_never_requested_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    
    if f"cum_nb_rejectedbps_per_{group_column}" in seqsum_df.columns:
        ax = plot_number_rejectedbps_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_number_sequenced_bps_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_number_ref_bps_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    
    ax = plot_read_duration_per_group(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_read_length_per_group_per_decisiontype(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    if "nb_ref_bps_full" in seqsum_df.columns:
        ax = plot_fraction_read_per_group_per_rejected(subsampled_seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(ax.figure)
    
    # require full seqsum_df
    fig = plot_number_channels_per_group_over_time(seqsum_df, save_dir=save_dir, group_column=group_column); logger.debug("Created 1 plot"); close_fig(fig)
    ax = plot_channel_occupation_fraction_over_time(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_channels_over_time(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(ax.figure)
    if "mux" in seqsum_df.columns and (seqsum_df["mux"].nunique() > 1):
        ax = plot_mux_over_time(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(ax.figure)
    fig, _ = plot_read_stats_by_channel_hists(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(fig)
    ax = plot_fraction_states_per_channel(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(ax.figure)
    ax = plot_read_end_reason_hist(seqsum_df, save_dir=save_dir); logger.debug("Created 1 plot"); close_fig(ax.figure)

def plot_coverage_per_group(cov_df, cov_thresholds=[1, 2, 3, 4, 5, 6], save_dir: Optional[Path]=None, group_column="group", close_figures=None):
    """Plot fraction covered per group for each coverage, then per coverage for each group"""
    if close_figures is None:
        close_figures = save_dir is not None
        
    def close_fig(fig):
        if close_figures:
            plt.close(fig)
    
    # fraction covered, one plot per coverage
    for cov in cov_thresholds:
        fig, ax = plt.subplots()
        sns.lineplot(cov_df, x="end_time", y=f"fraction_cov{cov}", hue=group_column, ax=ax)
        ax.set_xlabel("Read end time (s)")
        ax.set_ylabel("Fraction covered")
        ax.set_title(f"Coverage {cov}")
        if save_dir is not None:
            save_fig_and_pickle(fig, save_dir / f"fraction_cov{cov}_per_{group_column}.{FIGURE_EXT}")
        logger.debug("Created 1 plot"); close_fig(ax.figure)
        
    # fraction covered, one plot per group
    for (group, df) in cov_df.groupby(group_column, observed=True):
        fig, ax = plt.subplots()
        # reset_index to get index() column, cannot use non-unique end_time since per end_time, the coverage per chrom is recorded, and reads may also finish at exactly the same time
        df = pd.wide_to_long(df.reset_index(), stubnames="fraction_cov", i="index", j="cov").reset_index()
        df["cov"] = df["cov"].astype("category")
        sns.lineplot(df, x="end_time", y="fraction_cov", hue="cov", ax=ax)
        ax.set_xlabel("Read end time (s)")
        ax.set_ylabel("Fraction covered")
        ax.set_title(f"{group_column} {group}")
        if save_dir is not None:
            save_fig_and_pickle(fig, save_dir / f"fraction_covered_{group}.{FIGURE_EXT}")
        logger.debug("Created 1 plot"); close_fig(ax.figure)
            
def create_plots_for_seqsum(seqsum_df, nrows=None, group_to_units: Dict[str, List[Any]]=None, group_column=None, 
                            ref_genome_path=None, paf_file=None, cov_thresholds=[1, 2, 3, 4, 5, 6], cov_every=1, 
                            save_dir=None, close_figures=None):
    """
    Create plots for a sequencing summary file
    
    Args:
        seqsum_df: path to sequencing summary file, or dataframe
        nrows: only read the first nrows reads
        group_to_units: dictionary {group_name: units} where units form a subset of the unique values in group_column; if None, groups have size 1
        group_column: column in sequencing summary file to group by; if "all", use one group called "all"; if None, use GROUP_COLUMN
        
        ref_genome_path: path to reference genome; if None, don't plot coverage
        paf_file: path to PAF file to map reads to unit; if None, unit is the chromosome extracted from NanoSim read id
        cov_thresholds: coverage thresholds to plot
        cov_every: coverage is calculated every cov_every reads
        
        save_dir: directory to save plots to, if None, plots are not saved
        close_figures: close figures after saving, if None, close figures if save_dir is not None
        
    Returns:
        seqsum_df, cov_df
    """
    group_column = group_column or GROUP_COLUMN
        
    if save_dir is not None:
        save_dir.mkdir(exist_ok=True)
    
    if not isinstance(seqsum_df, pd.DataFrame):
        logger.debug(f"Reading {nrows if nrows is not None else 'all'} reads from sequencing summary file '{seqsum_df}'")
        seqsum_df_filename = seqsum_df
        seqsum_df = pd.read_csv(seqsum_df_filename, sep="\t", nrows=nrows)
        logger.debug(f"Done reading sequencing summary file '{seqsum_df_filename}'")
    
    if len(seqsum_df) == 0:
        logger.warning(f"Empty sequencing summary file '{seqsum_df}'")
        return seqsum_df, None
    
    logger.info(f"Sorting and cleaning seqsummary file of shape {seqsum_df.shape}")
    seqsum_df = sort_and_clean_seqsum_df(seqsum_df)
    logger.info(f"Adding previous gap duration to seqsummary")
    seqsum_df = add_previous_gap_duration(seqsum_df, seq_start_time=0)
    if group_column == "all":
        logger.info(f"Adding group column 'all'")
        seqsum_df[group_column] = "all"
        seqsum_df["nb_ref_bps"] = seqsum_df["sequence_length_template"]
    elif paf_file is not None:
        logger.info(f"Adding group column from PAF file '{paf_file}'")
        seqsum_df = add_group_and_reflen_from_paf(seqsum_df, paf_file=paf_file, group_column=group_column)
    elif group_column not in seqsum_df.columns:
        logger.info(f"Adding group column from NanoSim read id")
        add_group_and_reflen_from_nanosim_id(seqsum_df, group_column=group_column)
    
    chrom_column = group_column
    if group_to_units is not None:
        # create column to group by, e.g. several chromosomes in one group
        units_in_group = set.union(*[set(units) for units in group_to_units.values()])
        observed_units = set(seqsum_df[group_column].unique())
        if not units_in_group.issubset(observed_units):
            logger.warning(f"No reads were observed from the following groups: {units_in_group - observed_units}")
        other_group = observed_units - units_in_group
        if len(other_group) > 0:
            assert "other" not in group_to_units
            group_to_units["other"] = other_group
        logger.info(f"Plotting according to groups {group_to_units}")
        
        group_column = "group"
        assert group_column not in seqsum_df.columns, f"New column '{group_column}' already in sequencing summary df with columns {seqsum_df.columns}"
        unit_to_group = {unit: group for (group, units) in group_to_units.items() for unit in units}
        seqsum_df[group_column] = seqsum_df[chrom_column].map(unit_to_group)
        seqsum_df[group_column].fillna("other", inplace=True)
    logger.info("Adding extra columns for plotting")
    seqsum_df = seqsum_add_cols_for_plotting_selseq_performance(seqsum_df, group_column=group_column)
    
    logger.debug("Creating plots for seqsum...")
    plot_processed_seqsum(seqsum_df, group_column=group_column, save_dir=save_dir, close_figures=close_figures)
    logger.debug("Done creating plots for seqsum...")
    
    # create coverage-related plots
    if (ref_genome_path is not None) or (paf_file is not None):
        if paf_file is None:
            logger.debug(f"Computing coverage with respect to reference genome '{ref_genome_path}' for thresholds {','.join(map(str, cov_thresholds))} every {cov_every} reads")
            cov_tracker = NanoSimCoverageTracker.empty_from_ref_genomes(ref_genome_path, blocksize=1_000) #todo2: blocksize depending on chromosome size
        else:
            logger.debug(f"Computing coverage given PAF file '{paf_file}' for thresholds {','.join(map(str, cov_thresholds))} every {cov_every} reads")
            # extract lengths from the PAF file
            cov_tracker = PafCoverageTracker.empty_from_paf(paf_file=paf_file, blocksize=1_000)
            # cov_tracker = PafCoverageTracker.empty_from_ref_genomes(ref_genome_path, paf_file=paf_file, blocksize=1_000)
        cov_df = compute_coverage_per_group_df(seqsum_df, cov_tracker=cov_tracker, cov_thresholds=cov_thresholds, coverage_every=cov_every, group_to_chroms=group_to_units, chrom_column=chrom_column)
        logger.debug("Creating plots for coverage...")
        plot_coverage_per_group(cov_df, cov_thresholds=cov_thresholds, save_dir=save_dir, group_column="group", close_figures=close_figures)
        logger.debug("Done creating plots for coverage...")
    else:
        cov_df = None
    
    if save_dir is None:
        logger.debug("Displaying plots...")
        import signal; signal.signal(signal.SIGINT, signal.SIG_DFL) # to have plt respond to Ctrl+C, see https://gist.github.com/djwbrown/3e24bf4e0c5e9ee156a5?permalink_comment_id=4128738#gistcomment-4128738
        plt.show()
    else:
        logger.debug(f"Plots saved to '{save_dir.resolve()}'")
        
    return seqsum_df, cov_df
        
def main():
    """
    CLI entrypoint to create plots from a sequencing summary file
    """
    add_comprehensive_stream_handler_to_logger(logging.DEBUG)
    
    if is_test_mode():
        args = argparse.Namespace()
        args.ref_genome_path = "runs/enrich_usecase/chm13v2.0_normalized1000000firsttwo.fa.gz"
        args.seqsummary_filename = "runs/enrich_usecase/sequencing_summary.txt"
        # args.save_dir = None
        args.save_dir = Path("runs/enrich_usecase/figures")
        args.nrows = 4000
        args.cov_every = 10
        args.paf_file = None
    else:
        parser = argparse.ArgumentParser(description="Plotting script for sequencing summary file from the simulator")
        parser.add_argument("ref_genome_path", type=Path, help="Path to the reference genome")
        parser.add_argument("seqsummary_filename", type=Path, help="Path to the sequencing summary file")
        parser.add_argument("--nrows", type=int, default=None, help="Number of rows to read from the sequencing summary file")
        parser.add_argument("--save_dir", type=Path, default=None, help="Directory to save plots to")
        parser.add_argument("--cov_thresholds", type=str, default="1,2,3,4,5,6", help="Comma-separated list of target coverages to plot; set to '' if non-NanoSim reads")
        parser.add_argument("--targets", type=str, default=None, help="if provided, a comma-separated list, e.g. 'chr1,chr2'. Creates two groups on the plot, opposing targets to the rest")
        parser.add_argument("--cov_every", type=int, default=1, help="Compute coverage every x reads")
        parser.add_argument("--paf_file", type=Path, default=None, help="Path to the PAF file to use for coverage computation and to group by")
        
        args = parser.parse_args()
        print_args(args, logger=logger)
    
    cov_thresholds = list(map(int, args.cov_thresholds.split(",")))
    if args.targets is None:
        group_units = None
    else:
        group_units = {"targets": args.targets.split(",")}
    create_plots_for_seqsum(seqsum_df=args.seqsummary_filename, nrows=args.nrows, group_to_units=group_units, ref_genome_path=args.ref_genome_path, paf_file=args.paf_file, cov_thresholds=cov_thresholds, cov_every=args.cov_every, save_dir=args.save_dir)
    
    logger.debug("Done with plotting script")
    
if __name__ == "__main__":
    main()