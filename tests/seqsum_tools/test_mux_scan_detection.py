import numpy as np
import pandas as pd
from simreaduntil.seqsum_tools.mux_scan_detection import check_mux_scan_windows, find_mux_scans, remove_mux_scans

from simreaduntil.seqsum_tools.seqsum_preprocessing import sort_and_clean_seqsum_df


def test_mux_scan_removal(shared_datadir):
    sequencing_summary_file = shared_datadir / "zymo_short_seqsum.txt"
    
    df_read = pd.read_csv(sequencing_summary_file, sep="\t")
    seqsum_df = sort_and_clean_seqsum_df(df_read, min_gap_duration=None)#, min_gap_duration=0.05)

    mux_scan_windows = find_mux_scans(seqsum_df["start_time"].values, seqsum_df["end_time"].values, seqsum_df["mux"].values)
    
    assert mux_scan_windows.shape == (3, 2) # extracted seqsummary to have 3 mux scans

    check_mux_scan_windows(seqsum_df, mux_scan_windows)

    seqsum_df, mux_scan_boundaries = remove_mux_scans(seqsum_df, mux_scan_windows)

    assert all(np.diff(mux_scan_boundaries) > 0), "mux scan boundaries should be sorted and have positive length"

    # first value of each group is NA
    assert all(seqsum_df.groupby("channel", observed=True)["start_time"].diff().dropna() >= 0), "start_time should be increasing" 
    assert all(seqsum_df["start_time"] <= seqsum_df["end_time"]), "start_time should be <= end_time"
    assert all(seqsum_df.groupby("channel", observed=True)["end_time"].diff().dropna() >= 0), "start_time should be increasing" 
    assert all(seqsum_df.groupby("channel", observed=True)["nb_scans_before"].diff().dropna() >= 0), "mux scan index should be increasing"
