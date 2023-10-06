from contextlib import contextmanager
import logging
import numpy as np
import pandas as pd
import pytest
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, remove_comprehensive_stream_handler

from simreaduntil.shared_utils.nanosim_parsing import NanoSimId

@pytest.helpers.register
def get_random_seqsum_df():
    """
    Generate random seqsummary df by first generating read durations, then gaps between reads
    """
    random_state = np.random.default_rng(4)
        
    nb_reads = 1000
    channels = 1 + random_state.choice(6, size=nb_reads)
    read_durations = random_state.uniform(1.5, 10., size=nb_reads)
    read_durations = random_state.uniform(0.5, 1.5, size=nb_reads) + np.maximum(0, random_state.normal(5, 2.5))
    short_gap_durations = random_state.uniform(0.1, 0.5, size=(nb_reads+1)//2)
    long_gap_durations = random_state.uniform(15.1, 50.5, size=nb_reads - (nb_reads+1)//2)
    gap_durations = np.concatenate((short_gap_durations, long_gap_durations))
    random_state.shuffle(gap_durations) # in-place
    
    # first assign cumulative, then make each channel start at 0
    seqsum_df = pd.DataFrame({"channel": channels, "duration": read_durations})
    seqsum_df["start_time"] = (read_durations + gap_durations).cumsum()
    seqsum_df["start_time"] -= seqsum_df.groupby("channel", observed=True)["start_time"].transform("min")
    seqsum_df["start_time"] += random_state.uniform(0.05, 0.2, size=nb_reads) # to avoid 0
    seqsum_df["end_time"] = seqsum_df["start_time"] + seqsum_df["duration"]
    template_offset = random_state.uniform(0, 0.1*2/3, size=nb_reads)
    seqsum_df["template_start"] = seqsum_df["start_time"] + template_offset
    seqsum_df["template_duration"] = read_durations - template_offset
    seqsum_df["nb_ref_bps_full"] = np.rint((seqsum_df["duration"] - (seqsum_df["template_start"] - seqsum_df["start_time"])) * random_state.uniform(400, 500, size=nb_reads)).astype(np.int64)
    seqsum_df["passes_filtering"] = True
    
    # for plotting the seqsum, add rejections
    seqsum_df["chrom"] = random_state.choice(["chr1", "chr2"], size=nb_reads, p=[0.7, 0.3])
    
    # seqsum_df["end_reason"] = random_state.choice(["signal_positive", "data_service_unblock_mux_change"], size=nb_reads, p=[0.9, 0.1])
    n_chr1 = sum(seqsum_df["chrom"] == "chr1")
    seqsum_df.loc[seqsum_df["chrom"] == "chr1", "end_reason"] = random_state.choice(["signal_positive", "data_service_unblock_mux_change"], size=n_chr1, p=[0.2, 0.8])
    seqsum_df.loc[seqsum_df["chrom"] == "chr2", "end_reason"] = random_state.choice(["signal_positive", "data_service_unblock_mux_change"], size=nb_reads - n_chr1, p=[0.95, 0.05])
    
    seqsum_df["nb_ref_bps"] = seqsum_df["nb_ref_bps_full"]
    nb_ref_bps_full_when_unblocked = seqsum_df.loc[seqsum_df["end_reason"] == "data_service_unblock_mux_change", "nb_ref_bps_full"]
    seqsum_df.loc[seqsum_df["end_reason"] == "data_service_unblock_mux_change", "nb_ref_bps"] = np.rint(
        nb_ref_bps_full_when_unblocked * random_state.uniform(0.2, 0.9, size=len(nb_ref_bps_full_when_unblocked))
    ).astype(np.int64)
    # ensure at least one bp is rejected when the end reason is data_service_unblock_mux_change, already the case due to data generation process
    seqsum_df["sequence_length_template"] = np.rint(seqsum_df["nb_ref_bps"] * random_state.uniform(0.9, 1.1, size=nb_reads)).astype(np.int64)
    
    ref_lens_per_chrom = {"chr1": 1_000_000, "chr2": 1_000_000}
    def get_nanosim_id(chrom, read_nb, direction, ref_len):
        ref_pos = random_state.integers(0, ref_lens_per_chrom[chrom] - ref_len + 1)
        return str(NanoSimId(chrom=chrom, ref_pos=ref_pos, read_nb=read_nb, direction=direction, ref_len=ref_len))
    seqsum_df["read_id"] = list(map(get_nanosim_id, seqsum_df["chrom"], np.arange(nb_reads), random_state.choice(["F", "R"], size=nb_reads), seqsum_df["nb_ref_bps"]))
    
    seqsum_df.sort_values("end_time", inplace=True)

    return seqsum_df
