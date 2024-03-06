"""
Utility functions to help with the usecase
"""

import argparse
import functools
import logging
import os
from pathlib import Path
from textwrap import dedent
from typing import List, Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import bgzf
import toml
from simreaduntil.seqsum_tools.mux_scan_detection import remove_mux_scans_from_file, get_seqsum_filename_with_removed_mux_scans
from simreaduntil.seqsum_tools.seqsum_plotting import create_plots_for_seqsum
from simreaduntil.seqsum_tools.seqsum_preprocessing import add_previous_gap_duration, compute_median_pore_speed, sort_and_clean_seqsum_df
from simreaduntil.shared_utils.dna import get_sequence_names_in_fasta
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, setup_logger_simple

from simreaduntil.shared_utils.nanosim_parsing import case_convert_dna
from simreaduntil.shared_utils.plotting import ignore_tight_layout_warning
from simreaduntil.shared_utils.utils import delete_dir_if_exists, dill_dump, dill_load, force_eval_generator_function, num_lines_in_file, print_args, tqdm_with_name
from simreaduntil.shared_utils.dna import get_reverse_complement, get_random_DNA_seq
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId
from simreaduntil.shared_utils.utils import print_cmd_and_run
from simreaduntil.simulator.channel_stats import ChannelStats, plot_channel_stats
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.gap_sampling.gap_sampler_per_window_until_blocked import GapSamplerPerWindowUntilBlocked
from simreaduntil.simulator.gap_sampling.gap_sampling import RandomGapSampler
from simreaduntil.simulator.gap_sampling.inactive_active_gaps_replication import SingleChannelInactiveActiveReplicator
from simreaduntil.simulator.gap_sampling.rolling_window_gap_sampler import RollingWindowGapSamplerPerChannel
from simreaduntil.simulator.simulator import get_simulator_delay_over_time_df, plot_sim_actions, plot_simulator_delay_over_time
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.simulator.utils import set_package_log_level
from simreaduntil.usecase_helpers.readfish_plotting import get_chunk_mapping_time_over_time_df, get_chunk_wait_time_over_time_df, get_extra_basecall_delay_over_time_df, get_processing_time_per_read_over_time_df, get_throttle_over_time_df, plot_chunk_mapping_time, plot_chunk_waiting_time, plot_extra_basecalling_delay_per_iter, plot_readfish_processing_time, plot_throttle_over_time

logger = setup_logger_simple(__name__)
"""module logger"""

def random_nanosim_reads_gen(random_state=np.random.default_rng(2), length_range=(8, 17)):
    """
    Generate random reads with a valid NanoSim id; the reads will however not map to the reference genome
    
    Args:
        random_state: random state; if not provided, repeated calls to this function return different values
        length_range: range of read lengths
        
    Yields:
        (read_id, seq), where read_id is a valid NanoSim id
    """
    length_fcn = lambda: random_state.integers(*length_range)
    
    i = 0
    while True:
        i += 1
        
        chrom = random_state.choice(["chr1", "chr2"])
        direction = random_state.choice(["F", "R"])
        ref_len = length_fcn()
        assert ref_len > 0
        seq_len = ref_len
        ref_pos = random_state.integers(100, 200)
        
        read_id = NanoSimId(chrom=chrom, ref_pos=ref_pos, read_nb=i, direction=direction, ref_len=ref_len)
        
        yield (read_id, get_random_DNA_seq(seq_len, random_state=random_state))

# to load the FASTA file when the function is called rather than when the first read is requested (which may delay the simulation if an index has to be built first)
@force_eval_generator_function
def perfect_reads_gen(fasta_filename: Path, read_lens_range: Tuple[int], random_state=np.random.default_rng(1), nanosim_read_id=True):
    """
    Generate perfect reads that align to the reference genome
    
    Args:
        fasta_filename: filename
        read_lens_range: range of read lengths
        random_state: random state; if not provided, repeated calls to this function return different values
        nanosim_read_id: whether to return a valid NanoSim id as read_id
        
    Yields:
        (read_id, seq) where read_id is or is not a valid NanoSim id (read number starting at 1)
    """
    
    with pysam.FastaFile(fasta_filename) as fasta:
        ref_lengths = {ref: fasta.get_reference_length(ref) for ref in fasta.references}
        lens_str = ", ".join(f"{chrom}:{length:.2g}" for (chrom, length) in ref_lengths.items())
        logger.info(f"""Generating perfect reads from reference genome '{fasta.filename.decode()}' with total length {sum(ref_lengths.values()):.2g}, individual lengths are {lens_str}""")
        # ref_genomes = {ref: ref_fasta.fetch(ref) for ref in ref_fasta.references}
        # logger.info("Loaded reference genome")
        assert read_lens_range[1] - 1 <= min(ref_lengths.values()), f"read lengths {read_lens_range[1]-1} is long for reference lengths {ref_lengths.values()}"
        
        i = 0
        while True:
            i += 1
            # logger.debug(f"Generating {i}th read")
            chrom = random_state.choice(fasta.references)
            read_len = random_state.integers(*read_lens_range)
            start_pos = random_state.choice(ref_lengths[chrom] - read_len + 1) # position with respect to forward strand
            seq = fasta.fetch(chrom, start_pos, start_pos + read_len)
            # seq = ref_genomes[chrom][start_pos: start_pos + read_len] # may be inefficient depending on implementation
            direction = random_state.choice(["F", "R"])
            if direction == "R":
                seq = get_reverse_complement(seq)
            read_id = NanoSimId(chrom=chrom, ref_pos=start_pos, read_nb=f"nb{i}", direction=direction, ref_len=read_len) if nanosim_read_id else f"read{i}"
            # logger.debug(f"Generated {i}th read with id {read_id}")
            yield (str(read_id), seq)

def random_reads_gen(random_state=np.random.default_rng(2), length_range=(8, 17)):
    """
    Generator returning random reads with random ids
    
    Args:
        random_state: random state to use
        length_range: range of read lengths to generate
    """
    length_fcn = lambda: random_state.integers(*length_range)
    
    i = 0
    while True:
        i += 1
        yield (f"read{i}", get_random_DNA_seq(length_fcn(), random_state=random_state))

def get_fasta_open_method(filename, write: bool = False):
    """
    Get the open method for a file; gzipped files are opened with gzip.open
    
    Args:
        filename: path to file
        write: whether to open for writing
        
    Returns:
        open method
    """
    mode = "w" if write else "r"
    if str(filename).endswith(".gz"):
        logger.info(f"Detected .gz extension for file '{filename}'")
        # bgzip compression needed for pysam
        from Bio import bgzf
        return lambda: bgzf.open(filename, mode=f"{mode}t")
    else:
        return lambda: open(filename, mode=mode)
    
def normalize_fasta_gen(fasta_filename, seq_names=None, random_state=np.random.default_rng(2)):
    """
    Generator to normalize FASTA to upper case and replacing N letters
    
    Args:
        fasta_filename: path to fasta, can be gzipped
        seq_names: sequences to include; if None, all
        random_state: random state; if not provided, repeated calls to this function return different values
        
    Yields:
        (sequence name, SeqIO.SeqRecord)
    """
    
    # note: checked that both methods produce the same files using bash's cmp function
    # method = "pyfastx" # if seq_names is not None else "biopython", pyfastx broken, see my comment on https://github.com/lmdu/pyfastx/issues/63
    method = "biopython" # if seq_names is not None else "biopython"
    if method == "biopython":
        # takes 3m50s on human chr22 (including writing)
        with get_fasta_open_method(fasta_filename)() as f:
            # use SeqIO rather than pysam.FastaFile because pysam discards the additional chromosome name
            for (i, record) in enumerate(SeqIO.parse(f, "fasta")):
                if seq_names is None or (record.id in seq_names):
                    yield (record.id, SeqIO.SeqRecord(seq=Seq(case_convert_dna(record.seq, random_state=random_state)), id=record.id, description=f"{record.description} normalized"))
                    # print(record.format("fasta"))
    else:
        # takes 45s on human chr22 (including writing), probably faster because it can access the chromosome out-of-order using the offset unlike biopython which does not have a cache
        import pyfastx
        for record in pyfastx.Fasta(str(fasta_filename), full_name=True):
            description = record.description
            # allow either match to fasta header or just sequence name
            if seq_names is None or description in seq_names or record.name in seq_names:
                yield (description, SeqIO.SeqRecord(seq=Seq(case_convert_dna(record.seq, random_state=random_state)), id=record.name, description=f"{description} normalized"))

def normalize_fasta(in_fasta_path, out_fasta_path, seq_names=None, random_state=np.random.default_rng(2)):
    """
    Normalize a FASTA file to upper case and replacing N letters
    
    Args:
        in_fasta_path: path to input FASTA file
        out_fasta_path: path to output FASTA file
        seq_names: sequences to include; if None, all
        random_state: random state; if not provided, repeated calls to this function return different values
    """
    in_fasta_path = Path(in_fasta_path)
    out_fasta_path = Path(out_fasta_path)
    
    sequence_names_in = get_sequence_names_in_fasta(in_fasta_path)
    logger.info(f"Sequence names in: {sequence_names_in}" + (f", restricting to seqnames {seq_names}" if seq_names is not None else ""))
    if (seq_names is not None) and (not set(seq_names).issubset(sequence_names_in)):
        logger.warning(f"Requesting sequences not in sequence file: {set(seq_names).difference(sequence_names_in)}")
    
    # todo2: can be parallelized, profile first
    with get_fasta_open_method(out_fasta_path, write=True)() as f:
        SeqIO.write(tqdm_with_name(normalize_fasta_gen(in_fasta_path, seq_names=seq_names, random_state=random_state)), f, "fasta")
    logger.info(f"Written to '{out_fasta_path.resolve()}'")
    
    logger.info(f"Sequence names out: {get_sequence_names_in_fasta(out_fasta_path)}")
    
def normalize_fasta_cmd():
    """
    Command line script to normalize a FASTA file and write it back to a file
    """
    
    add_comprehensive_stream_handler_to_logger(None, level=logging.INFO)
    
    # argparse to parse in_fasta_path, out_fasta_path
    parser = argparse.ArgumentParser(description="Normalize a FASTA file to upper case and replacing N letters")
    parser.add_argument("in_fasta_path", type=Path, help="input FASTA file")
    parser.add_argument("out_fasta_path", type=Path, help="output FASTA file")
    parser.add_argument("--seq_names", type=str, help="restrict to a subset of sequences if provided, e.g. 'chr20,chr21'", default=None)
    args = parser.parse_args()
    print_args(args, logger=logger)
    in_fasta_path = args.in_fasta_path
    out_fasta_path = args.out_fasta_path
    seq_names = args.seq_names.split(",") if args.seq_names is not None else None
    
    normalize_fasta(in_fasta_path, out_fasta_path, seq_names=seq_names)
    
    logger.info("Done normalizing genome")

def get_gen_nanosim_reads_cmd(nanosim_dir, nanosim_model_prefix, ref_genome_path, reads_dir, n_reads_per_sim, n_procs=4, perfect=False, use_slurm=False):
    """
    Get command to generate NanoSim reads in a separate process, slurm also supported

    Args:
        nanosim_dir: path to NanoSim directory
        reads_dir: directory to save reads to
        n_reads_per_sim: number of reads to generate per simulation
        n_procs: number of processes per simulation
        perfect: whether to generate perfect reads
        use_slurm: whether to use slurm and run a bunch of simulations in parallel

    Returns:
        command as str to execute with print_cmd_and_run or os.system
    """

    reads_prefix = Path(reads_dir) / (("perfect_" if perfect else "") + "reads_seed$seed")

    if use_slurm:
        script_header = dedent(rf"""
            #!/usr/bin/bash
            #SBATCH --job-name=nanosim_reads
            #SBATCH --cpus-per-task={n_procs+1}
            #SBATCH --time=4:00:00
            #SBATCH --output="slurm_%x_%A_%a.out"
            #SBATCH --mem=4500M
            #SBATCH --array=1-20:{n_procs} # otherwise same seeds reused by NanoSim
            
            source ~/.bashrc
            
            seed=$SLURM_ARRAY_TASK_ID
            echo "Seed: $seed"
            
            set -e
    """).strip()
    else:
        script_header = dedent(rf"""
            #!/usr/bin/bash
            seed=1
    """).strip()
    nanosim_command = script_header + os.linesep + dedent(rf"""

        conda run -n nanosim python -c "import HTSeq; print(HTSeq.__version__)"
        
        # cd <correct_dir>
        conda run -n nanosim \
            python "{Path(nanosim_dir) / "src/simulator.py"}" genome \
            --model_prefix "{nanosim_model_prefix}" \
            --ref_g "{ref_genome_path}" \
            -dna_type linear \
            --output "{reads_prefix}" \
            --number {n_reads_per_sim} \
            --seed "$seed" \
            --strandness 0.5 \
            --basecaller guppy \
            --aligned_rate "100%" \
            --num_threads "{n_procs}" \
            {"--perfect" if perfect else ""} \
            --no_error_profile \
            --no_flanking
            #; exit
    """).strip()

    # print(nanosim_command, file=open(f"runs/jobs/generate_nanosim_reads.sh", "w"))
    # print_cmd_and_run(nanosim_command, dry=True)
    # print("You need to run this code (starting from the appropriate directory)!")

    # print_cmd_and_run(nanosim_command, dry=True)

    return nanosim_command

def get_cleaned_seqsum_filename(seqsum_filename):
    """Filename to remove sequencing summary with removed mux scans, in same directory"""
    seqsum_filename = Path(seqsum_filename)
    # don't add at end via Path.stem since it may have two extensions (.txt.gz), so this keeps the extension
    return seqsum_filename.parent / ("cleaned_" + seqsum_filename.name)

def remove_mux_scans_and_clean_if_inexistent(sequencing_summary_file):
    """
    Remove mux scans from seqsum file, clean it and save to file, if no such file exists
    
    It creates two intermediate files, one for removed mux scans, the other after cleaning.
    
    Args:
        sequencing_summary_file: seqsum file from which to remove mux scans
        
    Returns:
        sequencing summary with mux scans removed (data.frame)
    """
    seqsum_file_cleaned = get_cleaned_seqsum_filename(sequencing_summary_file)
    if seqsum_file_cleaned.exists():
        logger.info(f"Cleaned seqsum file '{seqsum_file_cleaned}' already exists, loading it")
        return pd.read_csv(seqsum_file_cleaned, sep="\t")
    
    seqsum_file_no_mux_scans = get_seqsum_filename_with_removed_mux_scans(sequencing_summary_file)
    if seqsum_file_no_mux_scans.exists():
        logger.info(f"Seqsum file '{seqsum_file_no_mux_scans}' with removed mux scans already exists, skipping mux scan removal")
        logger.debug(f"Reading in seqsum file '{seqsum_file_no_mux_scans}' for parameter extraction")
        seqsum_df = pd.read_csv(seqsum_file_no_mux_scans, sep="\t")
    else:
        logger.debug(f"Removing mux scans from seqsum file '{sequencing_summary_file}'")
        seqsum_df, mux_scan_boundaries = remove_mux_scans_from_file(sequencing_summary_file, save_filename=get_seqsum_filename_with_removed_mux_scans(sequencing_summary_file))
        
    logger.debug(f"Sorting and cleaning seqsum file")
    seqsum_df = sort_and_clean_seqsum_df(seqsum_df, min_gap_duration=1e-8)
    seqsum_df = add_previous_gap_duration(seqsum_df, seq_start_time=0)
    
    logger.debug(f"Writing cleaned seqsum file to '{seqsum_file_cleaned}'")
    seqsum_df.to_csv(seqsum_file_cleaned, sep="\t", index=False)
    logger.debug(f"Wrote cleaned seqsum file to '{seqsum_file_cleaned}'")
    
    return seqsum_df

def create_simparams_if_inexistent(sim_params_filename, seqsum_param_extr_file, n_channels: Optional[int], gap_sampler_from_seqsum_fcn):
    """
    Create simparams from sequencing summary and save to file, if no such file exists
    
    Args:
        sim_params_filename: path to sim params file, to load from or to save to if inexistent
        seqsum_param_extr_file: sequencing summary file to use for parameter extraction
        n_channels: number of channels to simulate; if None, inferred from seqsum
        gap_sampler_from_seqsum_fcn: function that takes a seqsum_df and returns a function that creates gap samplers
        
    Returns:
        whether a params file was created
    """
    if sim_params_filename.exists():
        logger.info(f"Sim params file '{sim_params_filename}' already exists, skipping sim params creation")
        return False
    else:
        seqsum_df = remove_mux_scans_and_clean_if_inexistent(seqsum_param_extr_file)
        if n_channels is None:
            n_channels = seqsum_df["channel"].nunique()
        
        logger.debug(f"Extracting sim params for the gap sampler with {n_channels} channels")
        gap_sampler_maker = gap_sampler_from_seqsum_fcn(seqsum_df)

        random_state = np.random.default_rng(1) # todo2: one random state, also in on_sim
        sim_params = SimParams(
            gap_samplers={f"ch{i+1}": gap_sampler_maker(random_state=random_state)[1] for i in range(n_channels)},
            bp_per_second=compute_median_pore_speed(seqsum_df), min_chunk_size=200, default_unblock_duration=0.1, seed=0,
        )

        logger.debug(f"Saving sim_params to file '{sim_params_filename}'")
        dill_dump(sim_params, sim_params_filename)
        
        return True
    
def get_gap_sampler_method(gap_sampler_type, n_channels_full=None):
    """
    Get function to create a gap sampler given string description
    
    Args:
        gap_sampler_type: string description of gap sampler type
        n_channels_full: number of channels of original run (the sequencing summary file does not contain completely broken channels)
    """
    logger.info(f"Getting gap sampler for type {gap_sampler_type}")
    if gap_sampler_type == "constant_gaps":
        return ConstantGapsUntilBlocked.from_seqsum_df
    elif gap_sampler_type == "sampler_per_window":  
        return GapSamplerPerWindowUntilBlocked.from_seqsum_df
    elif gap_sampler_type == "sampler_per_rolling_window_channel":
        return functools.partial(RollingWindowGapSamplerPerChannel.from_seqsum_df, n_channels_full=n_channels_full)
    elif gap_sampler_type == "replication":
        return SingleChannelInactiveActiveReplicator.from_seqsum_df
    elif gap_sampler_type == "random":
        return RandomGapSampler.from_seqsum_df
    else:
        assert gap_sampler_type is None, f"unknown gap sampler type '{gap_sampler_type}'"

def plot_sim_stats(run_dir, figure_dir):
    """
    Plot statistics from the simulator
    
    Args:
        run_dir: directory containing simulator outputs, e.g. action_results.csv, channel_stats.dill, log file
        figure_dir: directory to save figures to, must exist
    """
    assert figure_dir.exists()
    
    action_results_filename = run_dir / "action_results.csv"
    if action_results_filename.exists():
        logger.info(f"Loading simulator action results from '{action_results_filename}'")
        action_results_df = pd.read_csv(action_results_filename, sep="\t")
        logger.info(f"Loaded simulator action results from '{action_results_filename}'")
        if len(action_results_df) > 0:
            plot_sim_actions(action_results_df, save_dir=figure_dir)
        else:
            logger.warning("No actions were performed during the simulation")

    channel_stats_filename = run_dir / "channel_stats.dill"
    if channel_stats_filename.exists():
        logger.debug("Plotting channel stats")
        channel_stats: List[ChannelStats] = dill_load(channel_stats_filename)
        plot_channel_stats(channel_stats, save_dir=figure_dir)
    
def create_figures(seqsum, run_dir, figure_dir=None, delete_existing_figure_dir=False, **kwargs):
    """
    Create figures from a finished run
    
    Args:
        seqsum: sequencing summary file or dataframe to plot
        run_dir: run_dir containing outputs from the simulator
        figure_dir: where to store figures; if None, run_dir/figures; will be created if inexistent
        delete_existing_figure_dir: whether to delete directory if it exists
        **kwargs: to pass to create_plots_for_seqsum
    """
    run_dir = Path(run_dir)
    if figure_dir is None:
        figure_dir = run_dir / "figures"
    figure_dir = Path(figure_dir)
    if delete_existing_figure_dir:
        delete_dir_if_exists(figure_dir)
    figure_dir.mkdir(exist_ok=True)
    
    with set_package_log_level(logging.DEBUG):
        logger.debug(f"Visualizing the simulation results and saving to '{figure_dir.resolve()}'")
        ignore_tight_layout_warning()
        # warnings.filterwarnings("ignore", message="Only plotting full reads")
        nb_reads = len(seqsum) if isinstance(seqsum, pd.DataFrame) else num_lines_in_file(seqsum) - 1
        seqsum_df, cov_df = create_plots_for_seqsum(seqsum, cov_every=max(1, int(nb_reads/100)), save_dir=figure_dir, **kwargs)
        
        plot_sim_stats(run_dir, figure_dir)
        
        
def plot_log_file_metrics(log_filename, save_dir=None):
    """Parse log file and plot metrics"""
    logger.info(f"Plotting metrics from log file '{log_filename}' and saving to '{save_dir}'")
    
    df = get_simulator_delay_over_time_df(log_filename)
    fig = plot_simulator_delay_over_time(df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)
    
    proc_df = get_processing_time_per_read_over_time_df(log_filename)
    fig = plot_readfish_processing_time(proc_df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)

    throttle_df = get_throttle_over_time_df(log_filename)
    fig = plot_throttle_over_time(throttle_df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)

    basecall_delay_df = get_extra_basecall_delay_over_time_df(log_filename)
    fig = plot_extra_basecalling_delay_per_iter(basecall_delay_df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)
    
    chunk_waiting_time_df = get_chunk_wait_time_over_time_df(log_filename)
    fig = plot_chunk_waiting_time(chunk_waiting_time_df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)
    
    chunk_mapping_time_df = get_chunk_mapping_time_over_time_df(log_filename)
    fig = plot_chunk_mapping_time(chunk_mapping_time_df, save_dir=save_dir); logger.debug("Created 1 plot"); plt.close(fig)
    
def extract_errfile_from_condor_jobad(jobad_filename):
    """Extract the file where stderr is redirected to from a condor jobad file"""
    # find first match of "Err = "
    with open(jobad_filename, "r") as f:
        for line in f:
            if line.startswith("Err = "):
                return line[6:].strip()[1:-1] # temporary hack to get rid of quotes; could use htcondor's pythonbindings instead: python package classad
    return None

# def plot_condor_log_file_metrics(save_dir=None):
#     """Parse log filename from condor job add, then plot metrics from log"""
#     jobad_filename = os.environ.get("_CONDOR_JOB_AD", None)
#     if jobad_filename is not None:
#         # parse log filename from condor job ad, then process it
#         log_filename = extract_errfile_from_condor_jobad(jobad_filename)
#         if log_filename is not None:
#             logger.info(f"Plotting metrics from condor log file '{log_filename}'")
#             plot_log_file_metrics(log_filename, save_dir=save_dir)
#         else:
#             logger.warning(f"Did not find log file in condor job ad: {jobad_filename}")
#     else:
#         logger.warning("Did not find condor job ad environment variable '_CONDOR_JOB_AD', cannot plot metrics")