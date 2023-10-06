"""
Plot readfish stats to assess performance bottlenecks.
"""


import datetime
import itertools
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from simreaduntil.shared_utils.plotting import FIGURE_EXT
from simreaduntil.shared_utils.merge_axes import save_fig_and_pickle

from simreaduntil.shared_utils.plotting import make_tight_layout


def get_extra_basecall_delay_over_time_df(log_filename):
    """Plot basecall delay over time."""
    # extract lines like the following
    # 2023-08-17 08:21:52,311 - Total basecalling extra wait time: 2.73e-05s for 163400 basepairs --- readfish_wrappers.py:140 (basecall_minknow) INFO ##

    MARKER = "Total basecalling extra wait time: "
    def parse_line(line):
        # return time of log entry, extra wait time, number of basepairs
        log_time, remaining = line.split(" - ", maxsplit=1)
        # convert log_time to time
        log_time = datetime.datetime.strptime(log_time, "%Y-%m-%d %H:%M:%S,%f")
        remaining = remaining.split(MARKER)[1]
        return log_time, float(remaining.split("s for ")[0]), int(remaining.split("s for ")[1].split(" basepairs")[0])

    # line = "2023-08-17 08:21:52,311 - Total basecalling extra wait time: 2.73e-05s for 163400 basepairs --- readfish_wrappers.py:140 (basecall_minknow) INFO ##"
    # parse_line(line)
    
    with open(log_filename) as f:
        info = [parse_line(line) for line in f if MARKER in line]
        # info = list(itertools.islice((parse_line(line) for line in f if MARKER in line), 100))
    df = pd.DataFrame.from_records(info, columns=["time", "extra_wait_time", "nb_basepairs"])
    if len(df) > 0:
        df["time"] = (df["time"] - df["time"].iloc[0]).dt.total_seconds()
    df["avg_extra_wait_time"] = df["extra_wait_time"] / df["nb_basepairs"]
    df["avg_extra_wait_time"].fillna(0, inplace=True)
        
    return df

def plot_extra_basecalling_delay_per_iter(df, save_dir=None):
    """
    Plots the extra basecalling delay per iteration.
    
    One iteration is a called to get_basecalled_read_chunks() which logs the total time spent due to basecalling.
    The basecalling starts right when the function is called, so if the processing of the basecalled data takes longer,
    there is no extra delay due to basecalling.
    """
    df = df.sample(min(len(df), 200))

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 4))
    # to show all in one plot
    # see https://matplotlib.org/3.4.3/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
    # fig, ax = plt.subplots()
    # twin1 = ax.twinx()
    # twin2 = ax.twinx()
    # twin2.spines.right.set_position(("axes", 1.2))

    sns.lineplot(df, x="time", y="extra_wait_time", ax=ax1)
    ax1.set_xlabel("Time (of iteration) (s)")
    ax1.set_ylabel("Extra wait time (whole iteration) (s)")
    ax1.set_title("Extra waiting time (whole iteration)")
    sns.lineplot(df, x="time", y="avg_extra_wait_time", ax=ax2)
    ax2.set_xlabel("Time (of iteration) (s)")
    ax2.set_ylabel("Average waiting time per basepair (s)")
    ax2.set_title("Average extra waiting time")
    sns.lineplot(df, x="time", y="nb_basepairs", ax=ax3)
    ax3.set_xlabel("Time (of iteration) (s)")
    ax3.set_ylabel("Number of called basepairs at iteration")
    ax3.set_title("Number of called basepairs at iteration")

    fig.suptitle("Extra delay due to basecalling (delaying ReadFish)")

    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"extra_basecall_delay.{FIGURE_EXT}")
    
    return fig

def get_processing_time_per_read_over_time_df(log_filename):
    MARKER = "ReadFish processing speed: "
    def parse_line(line):
        # return time of log entry, number of reads, time elapsed
        log_time, remaining = line.split(" - ", maxsplit=1)
        # convert log_time to time
        log_time = datetime.datetime.strptime(log_time, "%Y-%m-%d %H:%M:%S,%f")
        remaining = remaining.split(" --- ", maxsplit=1)[0]
        remaining = remaining.split(MARKER, maxsplit=1)[1]
        num_reads, time_elapsed = remaining.split("/")
        return log_time, int(num_reads[:-1]), float(time_elapsed[:-1])

    # line = "2023-08-17 08:21:38,824 - 51R/0.35265s --- ru_gen.py:388 (simple_analysis) INFO ##"
    # parse_line(line)
    
    with open(log_filename) as f:
        info = [parse_line(line) for line in f if MARKER in line]
        # info = list(itertools.islice((parse_line(line) for line in f if MARKER in line), 100))
    df = pd.DataFrame.from_records(info, columns=["time", "nb_reads", "elapsed_time"])
    if len(df) > 0:
        df["time"] = (df["time"] - df["time"].iloc[0]).dt.total_seconds()
    df["elapsed_time_per_read"] = df["elapsed_time"] / df["nb_reads"]
    df["elapsed_time_per_read"].fillna(0, inplace=True)
    
    return df

def plot_readfish_processing_time(df, save_dir=None):
    """Plot ReadFish processing time per read over time"""

    df = df.sample(min(len(df), 200))

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 4))
    # to show all in one plot
    # see https://matplotlib.org/3.4.3/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
    # fig, ax = plt.subplots()
    # twin1 = ax.twinx()
    # twin2 = ax.twinx()
    # twin2.spines.right.set_position(("axes", 1.2))

    sns.lineplot(df, x="time", y="elapsed_time", ax=ax1)
    ax1.set_xlabel("Time (of iteration) (s)")
    ax1.set_ylabel("Elapsed time (whole iteration) (s)")
    sns.lineplot(df, x="time", y="elapsed_time_per_read", ax=ax2)
    ax2.set_xlabel("Time (of iteration) (s)")
    ax2.set_ylabel("Average elapsed time per read (s)")
    sns.lineplot(df, x="time", y="nb_reads", ax=ax3)
    ax3.set_xlabel("Time (of iteration) (s)")
    ax3.set_ylabel("Number of reads at iteration")

    fig.suptitle("ReadFish processing time")

    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"readfish_stats_per_iter.{FIGURE_EXT}")

    return fig

def get_throttle_over_time_df(log_filename):
    MARKER = "ReadFish throttle: "
    def parse_line(line):
        # return time of log entry, throttle
        log_time, remaining = line.split(" - ", maxsplit=1)
        # convert log_time to time
        log_time = datetime.datetime.strptime(log_time, "%Y-%m-%d %H:%M:%S,%f")
        remaining = remaining.split(" --- ", maxsplit=1)[0]
        remaining = remaining.split(MARKER)[1]
        return log_time, float(remaining[:-1])
    
    # line = "2023-08-17 08:18:47,800 - Throttle: 0.00052s --- ru_gen.py:391 (simple_analysis) INFO ##"
    # parse_line(line)
    
    with open(log_filename) as f:
        info = [parse_line(line) for line in f if MARKER in line]
        # info = list(itertools.islice((parse_line(line) for line in f if MARKER in line), 100))
    df = pd.DataFrame.from_records(info, columns=["time", "throttle"])
    if len(df) > 0:
        df["time"] = (df["time"] - df["time"].iloc[0]).dt.total_seconds()
    
    return df

def plot_throttle_over_time(df, save_dir=None):
    """Plot ReadFish throttle over time"""

    df = df.sample(min(len(df), 200))
    df.sort_values("time", inplace=True)

    fig, ax = plt.subplots()
    ax.plot(df["time"], df["throttle"])
    sns.lineplot(df, x="time", y="throttle", ax=ax)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Throttle")
    ax.set_title("Throttle over time")
    
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"readfish_throttle.{FIGURE_EXT}")
    
    return fig


if __name__ == "__main__":
    log_filename = "/Volumes/mmordig/joblogs/job-13931078-0.err"

    proc_df = get_processing_time_per_read_over_time_df(log_filename)
    plot_readfish_processing_time(proc_df)

    throttle_df = get_throttle_over_time_df(log_filename)
    plot_throttle_over_time(throttle_df)

    basecall_delay_df = get_extra_basecall_delay_over_time_df(log_filename)
    plot_extra_basecalling_delay_per_iter(basecall_delay_df)
    