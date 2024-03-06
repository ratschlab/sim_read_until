"""
Combine the simulator with ReadFish (and optionally NanoSim)
"""

import argparse
from contextlib import contextmanager
import logging
from pathlib import Path
import signal
import time
import numpy as np
import pysam
from simreaduntil.shared_utils.debugging_helpers import is_test_mode
from simreaduntil.shared_utils.dna import get_ref_lengths

from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, setup_logger_simple
from simreaduntil.shared_utils.timing import cur_ns_time
from simreaduntil.shared_utils.utils import delete_dir_if_exists, dill_dump, dill_load, print_args, set_signal_handler
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.readpool import ReadPoolFromFile, ReadPoolFromIterable, ThreadedReadPoolWrapper
from simreaduntil.simulator.readswriter import ArrayReadsWriter, CompoundReadsWriter, RotatingFileReadsWriter, SingleFileReadsWriter, ThreadedReadsWriterWrapper
from simreaduntil.simulator.simfasta_to_seqsum import SequencingSummaryWriter, convert_simfasta_dir_to_seqsum
from simreaduntil.simulator.simulator import ONTSimulator, convert_action_results_to_df, run_periodic_mux_scan_thread, write_simulator_stats, stop_simulation_after_time_thread
from simreaduntil.simulator.simulator_client import DeviceGRPCClient
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.simulator.simulator_server import launchable_device_grpc_server, manage_grpc_server
from simreaduntil.simulator.utils import set_package_log_level

# export PYTHONPATH=".:$PYTHONPATH" to make these imports work
from simreaduntil.usecase_helpers.readfish_wrappers import BASECALLING_TIME_PER_BP, DummyBasecaller, NanoSimMapper, ReadUntilClientWrapper, get_flowcell_array_replacement, replace_ru_mapper, send_message_to_logger
from simreaduntil.usecase_helpers.utils import perfect_reads_gen

# logging setup, do this before readfish imports so that their NOTSET inherits from this root logger
logger = setup_logger_simple(__name__)
"""module logger"""

# uncomment this for IntelliSense to work in VSCode, e.g. for targeted_seq, see below ru imports as well
# # ReadFish assumes it is installed as ru package (through imports in the package, so we can't import it like this)
# import external.ont_readfish.ru.utils
# from external.ont_readfish.ru.unblock_all import simple_analysis as unblock_all
# from external.ont_readfish.ru.ru_gen import simple_analysis as targeted_seq, get_chunk_logger, get_paf_logger, get_run_info

# ru = ReadFish
import ru.utils
ru.utils.send_message = send_message_to_logger # make sure that ru.utils was not imported yet by another module in ru
ru.utils.get_flowcell_array = get_flowcell_array_replacement
from ru.unblock_all import simple_analysis as unblock_all # imports send_message, so should do this at the beginning
from ru.ru_gen import simple_analysis as targeted_seq, get_chunk_logger, get_paf_logger, get_run_info # imports send_message, so should do this at the beginning

def compute_nonselective_coverage(ref_genome_path, reads_file):
    """Compute coverage if all reads are played back without any selective sequencing, i.e. full reads"""
    ref_length = sum(get_ref_lengths(ref_genome_path).values())
    reads_file = Path(reads_file)
    if reads_file.is_dir():
        total_reads_length =  sum(sum(get_ref_lengths(x).values()) for x in reads_file.glob("*.fasta"))
    else:
        total_reads_length =  sum(get_ref_lengths(reads_file).values())
    return total_reads_length / ref_length

def get_reads_writer(run_dir: Path, rotating_writeout: bool):
    # reads_writer = ArrayReadsWriter() # for debugging mostly
    
    mk_run_dir = run_dir / "reads"
    delete_dir_if_exists(mk_run_dir)
    mk_run_dir.mkdir()

    if rotating_writeout:
        reads_writer = RotatingFileReadsWriter(mk_run_dir, "reads_", max_reads_per_file=4000)
    else:
        reads_writer = SingleFileReadsWriter(open(mk_run_dir / "reads.fasta", "w"))
    
    seqsum_writer = SequencingSummaryWriter(open(run_dir / "live_sequencing_summary.txt", "w"))
    reads_writer = CompoundReadsWriter([reads_writer, seqsum_writer])
    reads_writer = ThreadedReadsWriterWrapper(reads_writer)
    return reads_writer

def get_sim_params(sim_params_file, n_channels) -> SimParams:
    if sim_params_file is None:
        # take realistic params, otherwise ReadFish mapper (minimap2) will still not map after 12 (small) chunks
        sim_params = SimParams(
            gap_samplers={f"ch{i+1}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0.1) for i in range(n_channels)},
            bp_per_second=450, min_chunk_size=200, default_unblock_duration=0.1, seed=0,
        )
    else:
        logger.info(f"Loading simparams from '{sim_params_file}'")
        sim_params: SimParams = dill_load(sim_params_file)
        if n_channels != sim_params.n_channels:
            logger.warning(f"Using sim_params.n_channels={sim_params.n_channels} instead of {n_channels} because it was saved in the sim_params_file")
        
    assert sorted(list(sim_params.gap_samplers.keys())) == sorted([f"ch{i+1}" for i in range(n_channels)]) # assumed by downstream plotting scripts
    
    return sim_params
    
def get_read_pool(reads_file_type, reads_file, ref_genome_path, n_channels, reads_len_range=None):
    """Get read pool either from reads_file or perfect reads from ref_fasta"""
        
    if reads_file_type == "generate":
        # read_pool = ReadPoolFromIterable(random_nanosim_reads_gen(random_state=np.random.default_rng(3), length_range=(10, 50)))
        logger.info(f"Generating perfect reads without NanoSim using ref genome '{ref_genome_path}'")
        assert reads_len_range is not None
        read_pool = ReadPoolFromIterable(perfect_reads_gen(ref_genome_path, read_lens_range=reads_len_range, random_state=np.random.default_rng(1)))
    elif reads_file_type == "fasta":
        logger.info("Reading in FASTA reads. pysam index creation may take some time.")
        # read_pool = ReadPoolFromFile(reads_file=reads_file)
        # read_pool = ReadPoolFromFile(reads_file_or_dir=reads_file)
        read_pool = ReadPoolFromFile(reads_file_or_dir=reads_file, shuffle_rand_state=np.random.default_rng(3)) # use a different rng since not thread-safe!
    else:
        logger.info("Creating slow5 read pool")
        raise "slow5 currently unavailable"
        # read_pool = Slow5ReadPool(reads_file, read_buffer=2*n_channels) # todo: remove read_buffer
        
    read_pool = ThreadedReadPoolWrapper(read_pool, queue_size=2*n_channels)
    return read_pool

@contextmanager
def wrap_simulator_in_grpc(simulator, use_grpc):
    """Launch a gRPC server proxying the simulator, if use_grpc is set"""
    if use_grpc:
        port, server, unique_id = launchable_device_grpc_server(simulator)
        assert port != 0
        with manage_grpc_server(server):
            with DeviceGRPCClient(port) as client:
                yield client
    else:
        # don't do anything
        yield simulator

"""
Auto-detect the reads file type

Args:
    reads_file (Path): path to the reads file or directory, possibly None
    
Returns:
    If reads file is None, "generate:
    Else If reads file ends with fasta or is a directory containing fastas, "fasta"
    Else If reads file ends with slow5/blow5 or is a directory containing slow5/blow5, "slow5"
    Else raise ValueError
"""
def get_reads_file_type(reads_file: Path):
    if reads_file is None:
        return "generate" # generate reads
    elif ReadPoolFromFile.can_handle(reads_file):
        return "fasta"
    # elif Slow5ReadPool.can_handle(reads_file):
    #     return "slow5"
    else:
        raise ValueError(f"Cannot determine reads file type for file '{reads_file}'")
    
    # if reads_file is None:
    #     return "generate" # generate reads
    # elif reads_file.endswith(".fasta"):
    #     return "fasta"
    # else:
    #     return "slow5"
    
def main(toml_file):
    """
    Run ReadFish with the simulator
    
    It outputs a paf log file and a chunk log file in the run_dir.
    The paf log file contains a mapping each time a read chunk is mapped.
    The chunk log file contains information about the chunk decisions.
    
    Args:
        toml_file (Path): toml file from which to load parameters
        
    Returns:
        path to the sequencing summary file, None if not written
    """
    
    ####################################################
    #################### PARAMETERS ####################
    ####################################################

    import toml
    config = toml.load(toml_file)

    run_dir = Path(config["run_dir"])
    acceleration_factor = float(config["acceleration_factor"]) if "acceleration_factor" in config else 1
    realtime_run_duration = float(config["run_duration"]) / acceleration_factor
    
    reads_file = config.get("reads_file", None)
    reads_len_range = config.get("reads_len_range", None)
    if reads_len_range is not None:
        assert isinstance(reads_len_range, list)
    ref_genome_path = config.get("ref_genome_path", None)
    sim_params_file = config.get("sim_params_file", None)
    rotating_writeout = config.get("rotating_writeout", True)
    use_grpc = config.get("use_grpc", False)
    if "mux_scan_period" in config:
        mux_scan_period = float(config["mux_scan_period"])
        mux_scan_duration = float(config["mux_scan_duration"])
    else:
        mux_scan_period = None
        mux_scan_duration = None
    write_seqsum_file = config.get("write_seqsum_file", True)

    # readfish params
    readfish_config_file = Path(config["readfish_config_file"]) if "readfish_config_file" in config else None
    readfish_method = config.get("readfish_method", "targeted_seq")
    assert readfish_method in ["control", "unblock_all", "targeted_seq"]
    logger.info(f"Running ReadFish with method '{readfish_method}'")
    
    ####################################################
    #################### SET UP RUN ####################
    ####################################################

    reads_file_type = config.get("reads_file_type", get_reads_file_type(reads_file))
    logger.info(f"Using reads file type '{reads_file_type}'")
    if (ref_genome_path is not None) and reads_file_type == "fasta":
        logger.info(f"(Without selseq,) You will reach average coverage {compute_nonselective_coverage(ref_genome_path, reads_file=reads_file)}")

    if run_dir.exists():
        logger.warning(f"Run dir '{run_dir}' already exists")
    run_dir.mkdir(exist_ok=True)
    sim_params = get_sim_params(sim_params_file=sim_params_file, n_channels=int(config["n_channels"]))
    reads_writer = get_reads_writer(run_dir, rotating_writeout=rotating_writeout)
    read_pool = get_read_pool(reads_file_type=reads_file_type, reads_file=reads_file, ref_genome_path=ref_genome_path, n_channels=sim_params.n_channels, reads_len_range=reads_len_range)
    # if isinstance(read_pool, Slow5ReadPool):
    #     logger.info("Detected slow5 reads, so adapting chunk size and bp_per_second")
    #     sim_params.bp_per_second = int(sim_params.bp_per_second * 4000/450) # todo: rename bp_per_second
    #     sim_params.min_chunk_size = int(sim_params.min_chunk_size * 4000/450)

    # When the simulation is accelerated, we decrease the batch size so that ReadFish can send out actions faster. 
    # Higher acceleration factor means that reads will finish faster.
    # The throttle is decreased because there are less reads per batch and the acceleration is going faster.
    original_batch_size = int(config.get("readfish_batch_size", sim_params.n_channels))
    # readfish_batch_size = max(1, round(original_batch_size / acceleration_factor))
    readfish_batch_size = original_batch_size
    readfish_throttle = float(config.get("readfish_throttle", 0.1)) * (readfish_batch_size / original_batch_size) / acceleration_factor
    logger.info(f"Running ReadFish with batch size {readfish_batch_size} and throttle {readfish_throttle}s")
    
    simulator : ONTSimulator = ONTSimulator(
        read_pool=read_pool, 
        reads_writer=reads_writer,
        sim_params=sim_params,
        output_dir=run_dir,
    )
    
    def start_sim():
        # start sim and mux scan thread
        logger.info(f"Starting the simulation with {simulator.n_channels} channels")
        simulator.start(acceleration_factor=acceleration_factor)
        logger.info(f"Started the simulation")
        if mux_scan_period is not None:
            mux_scan_thread = run_periodic_mux_scan_thread(simulator, period=mux_scan_period, scan_duration=mux_scan_duration, acceleration_factor=acceleration_factor)
            mux_scan_thread.start()
        else:
            mux_scan_thread = None
        return mux_scan_thread
    
    orig_simulator = simulator
    with wrap_simulator_in_grpc(orig_simulator, use_grpc=use_grpc) as simulator:
        # use a thread to stop it because simulator.stop() may never be able to acquire the mutex in the signal handler because the mutex may still be held by simulator.forward() or similar when the signal handler is running
        with set_signal_handler(signal_type=signal.SIGINT, handler=lambda *args, **kwargs: stop_simulation_after_time_thread(simulator, t=0).start()):
            if readfish_method == "control":
                mux_scan_thread = start_sim()
                
                stop_thread = stop_simulation_after_time_thread(simulator, t=realtime_run_duration)
                stop_thread.start()
                
            elif readfish_method == "unblock_all":
                mux_scan_thread = start_sim()
            
                stop_thread = stop_simulation_after_time_thread(simulator, t=realtime_run_duration)
                stop_thread.start()
                
                unblock_all(ReadUntilClientWrapper(simulator), duration=realtime_run_duration, batch_size=readfish_batch_size, throttle=readfish_throttle, unblock_duration=0.1)
                
            elif readfish_method == "targeted_seq":
                chunk_logger = get_chunk_logger(run_dir / "chunk_log.txt")
                paf_logger = get_paf_logger(run_dir / "mapping.paf")
                live_toml = Path("{}_live".format(readfish_config_file))
                dummy_caller = DummyBasecaller(time_per_bp=BASECALLING_TIME_PER_BP / acceleration_factor)
                run_info, conditions, reference, caller_kwargs = get_run_info(
                    readfish_config_file, num_channels=simulator.n_channels,
                )
                
                use_fake_mapper = (reference == "fake_mapper")
                with replace_ru_mapper(replace=use_fake_mapper) as CustomMapper:
                    # load Minimap2 index
                    if not use_fake_mapper: logger.info(f"Initializing minimap2 mapper from file '{Path(reference).resolve()}'")
                    t0 = cur_ns_time()
                    mapper = CustomMapper(reference)
                    if not use_fake_mapper: logger.info(f"Mapper initialized after {cur_ns_time() - t0}s")
                
                    # only start once the custom mapper is initialized (as in ReadFish ru_gen.py)
                    mux_scan_thread = start_sim()
                    
                    stop_thread = stop_simulation_after_time_thread(simulator, t=realtime_run_duration)
                    stop_thread.start()
                    
                    targeted_seq(ReadUntilClientWrapper(simulator), batch_size=readfish_batch_size, throttle=readfish_throttle, unblock_duration=0.1, 
                                chunk_logger=chunk_logger, paf_logger=paf_logger, live_toml_path=live_toml,
                                flowcell_size=simulator.n_channels,
                                dry_run=False, run_info=run_info, conditions=conditions, mapper=mapper, caller=dummy_caller)
                    
            stop_thread.raise_if_error()
            while simulator.is_running:
                # stop call() returns False early, when another stop call is already running
                time.sleep(0.5)
                    
            if mux_scan_thread is not None:
                mux_scan_thread.raise_if_error()
    simulator = orig_simulator # restore
    # go outside grpc connection
    write_simulator_stats([simulator], output_dir=simulator.mk_run_dir)
    
    read_pool.finish()
    
    seqsum_filename = None
    if isinstance(reads_writer, ArrayReadsWriter):
        # reads_writer.reads
        logger.info(reads_writer.extended_repr())
    else:
        # logger.info(reads_writer)
        
        # todo: seqsummary file written on the fly, remove write_seqsum_file arg,
        if write_seqsum_file:
            seqsum_filename = run_dir / "sequencing_summary.txt"
            logger.info(f"Writing sequencing summary file '{seqsum_filename}'")
            convert_simfasta_dir_to_seqsum(run_dir / "reads", seqsummary_filename=seqsum_filename)
            logger.info("Wrote sequencing summary file")

    logger.info("Done with simulation")
    
    return seqsum_filename

if __name__ == "__main__":
    log_level = logging.INFO
    logging.getLogger(__name__).setLevel(log_level)
    add_comprehensive_stream_handler_to_logger(None, level=log_level)
    with set_package_log_level(log_level):
        if is_test_mode():
            import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/dummy_replication").expanduser())
            main(Path("usecases/configs/readfish_with_sim_config.toml"))
        else:
            parser = argparse.ArgumentParser(description="Run ReadFish with the ONTSimulator")
            parser.add_argument("--toml", type=Path, help="toml file from which to load parameters", default=Path("configs/config.toml"))
            args = parser.parse_args()
            print_args(args, logger=logger)
            
            main(args.toml)
            