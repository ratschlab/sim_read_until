import pandas as pd
import pytest
from pytest import approx
from pandas.testing import assert_frame_equal

from simreaduntil.seqsum_tools.seqsum_preprocessing import add_previous_gap_duration, ensure_min_gaps_between_reads_single, get_gaps_single_channel

def test_ensure_min_gaps_between_reads_single():
    min_gap = 0.3
    df_single = pd.DataFrame([
        (2, 3), (4, 6.5), (7, 9), (20, 5)
    ], columns=["start_time", "duration"])
    df_single["template_start"] = df_single["start_time"]
    df_single["channel"] = 1

    expected_df = pd.DataFrame([
        (2, 3), (5.3, 6.5), (12.1, 9), (25.1, 5)
    ], columns=["start_time", "duration"])
    expected_df["template_start"] = expected_df["start_time"]
    expected_df["channel"] = 1
    
    assert_frame_equal(ensure_min_gaps_between_reads_single(df_single, min_gap_duration=min_gap), expected_df, check_like=True)
    
def test_add_previous_gap_duration():
    df = pd.DataFrame([
        [1, 3.0, 10.0],
        [1, 11.1, 17.3],
        #
        [2, 2.0, 3.0],
        [2, 6.0, 9.0],
    ], columns=["channel", "start_time", "end_time"])
    expected_df = pd.DataFrame([
        [1, 3.0, 10.0, 2.0],
        [1, 11.1, 17.3, 1.1],
        #
        [2, 2.0, 3.0, 1.0],
        [2, 6.0, 9.0, 3.0],
    ], columns=["channel", "start_time", "end_time", "prev_gap_duration"])

    assert_frame_equal(
        add_previous_gap_duration(df, seq_start_time=1.0),
        expected_df,
        check_like=True
    )

    with pytest.raises(AssertionError):
        # not sorted
        add_previous_gap_duration(df.iloc[[1, 0]], seq_start_time=1.0)
        
def test_get_gaps_single_channel():
    df_single = pd.DataFrame([
        (6, 9), (9.5, 10.5), (12.5, 13.5)
    ], columns=["start_time", "end_time"])
    df_single["channel"] = 1
    assert get_gaps_single_channel(df_single, seq_start_time=2) == (approx([4, 0.5, 2]), approx([2, 9, 10.5]))