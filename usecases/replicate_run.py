"""
Extract parameters from a sequencing summary, run it (without selective sequencing) and compare the statistics to the original run.

The simulation is replicated channelwise. It ignores the following config values: ref_genome_path, rotating, mux_scan*, use_grpc, readfish*.

cd runs/run_replication
source ~/ont_project_all/ont_project_venv/bin/activate
export PATH=~/ont_project_all/tools/bin:$PATH && which minimap2

mkdir random && cd random
ln -s ../data .
ln -s ../configs/random configs
python ~/ont_project_all/ont_project/usecases/replicate_run.py

cd ..; rm -rf sampler_per_rolling_window_channel
mkdir sampler_per_rolling_window_channel && cd sampler_per_rolling_window_channel
ln -s ../data .
ln -s ../configs/sampler_per_rolling_window_channel configs
python ~/ont_project_all/ont_project/usecases/replicate_run.py

"""

import logging
from pathlib import Path
import pandas as pd

import toml
from simreaduntil.seqsum_tools.seqsum_plotting import create_plots_for_seqsum
from simreaduntil.shared_utils.debugging_helpers import is_test_mode
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, print_logging_levels, setup_logger_simple
from simreaduntil.shared_utils.plotting import filter_seaborn_warnings, ignore_tight_layout_warning
from simreaduntil.shared_utils.utils import delete_dir_if_exists, dill_dump, dill_load, num_lines_in_file, subset_dict
from simreaduntil.simulator.gap_sampling.inactive_active_gaps_replication import get_read_durations_per_channel
from simreaduntil.simulator.simfasta_to_seqsum import convert_simfasta_dir_to_seqsum, convert_simfasta_to_seqsum
from simreaduntil.simulator.simulator import assign_read_durations_to_channels, run_simulator_from_sampler_per_channel, run_simulator_from_sampler_per_channel_parallel, write_simulator_stats
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.simulator.utils import set_package_log_level
from simreaduntil.usecase_helpers.utils import get_cleaned_seqsum_filename, create_figures, create_simparams_if_inexistent, get_gap_sampler_method, remove_mux_scans_and_clean_if_inexistent

logger = setup_logger_simple(__name__)

add_comprehensive_stream_handler_to_logger(None)
logging.getLogger(__name__).setLevel(logging.DEBUG)
logging.getLogger("simreaduntil").setLevel(logging.DEBUG)
# logging.getLogger().setLevel(logging.DEBUG) # warnings from everywhere, not desired

filter_seaborn_warnings()
print_logging_levels()

if is_test_mode():
    # todo
    # import os; os.chdir(Path("~/ont_project_all/ont_project/runs/run_replication/dummy_random").expanduser())
    import os; os.chdir(Path("~/ont_project_all/ont_project/runs/run_replication/sampler_per_rolling_window_channel").expanduser())
sim_config_file = Path("configs/config.toml")
ask_dir_deletion = False

sim_config = toml.load(sim_config_file)
run_dir = Path(sim_config["run_dir"])
# TeeStdouterr(run_dir / "stdouterr.txt").redirect()
logger.debug(f"Read in simulation config file '{sim_config_file}'")

sim_params_filename = Path(sim_config["sim_params_file"]) if "sim_params_file" in sim_config else None
seqsum_param_extr_file = Path(sim_config["seqsum_param_extr_file"]) if "seqsum_param_extr_file" in sim_config else None
gap_sampler_type = sim_config.get("gap_sampler_type", None)
gap_sampler_method = get_gap_sampler_method(gap_sampler_type, n_channels_full=sim_config["n_channels_full"])

logger.info("#"*80)
logger.info(f"Loaded config file '{sim_config_file}' with content:\n{sim_config_file.read_text()}")
logger.info("#"*80)

def get_read_durations_filename(seqsum_filename):
    """Filename to save read durations extracted from sequencing summary file"""
    seqsum_filename = Path(seqsum_filename)
    # don't add at end via Path.stem since it may have two extensions (.txt.gz), so this keeps the extension
    return seqsum_filename.parent / ("read_durations_" + seqsum_filename.name)

def extract_read_durations_if_inexistent(seqsum_filename):
    read_durations_filename = get_read_durations_filename(seqsum_filename)
    if read_durations_filename.exists():
        logger.debug(f"Read durations file '{read_durations_filename}' already exists, loading it")
        return dill_load(read_durations_filename)
    else:
        logger.debug(f"Extracting read durations from sequencing summary file '{seqsum_filename}'")
        seqsum_df = pd.read_csv(seqsum_filename, sep="\t") # same file as used for simparams
        read_durations_per_channel = get_read_durations_per_channel(seqsum_df)
        dill_dump(read_durations_per_channel, read_durations_filename)
        logger.debug(f"Saved read durations to '{read_durations_filename}'")
    return read_durations_per_channel

def run_simulator(seqsum_filename):
    create_simparams_if_inexistent(sim_params_filename, seqsum_param_extr_file, sim_config.get("n_channels", None), gap_sampler_method)
    sim_params: SimParams = dill_load(sim_params_filename)
    logger.debug(f"Loaded gap sampler for {sim_params.n_channels} channels")

    read_durations_per_channel = extract_read_durations_if_inexistent(get_cleaned_seqsum_filename(seqsum_param_extr_file))
    # possibly restrict if simulating less channels
    read_durations_per_channel = assign_read_durations_to_channels(read_durations_per_channel.values(), channel_names=sim_params.gap_samplers.keys())

    delete_dir_if_exists(run_dir, ask=ask_dir_deletion)
    reads_dir = run_dir / "reads"
    reads_dir.mkdir(parents=True)

    logger.debug(f"#################################################################")
    logger.debug(f"#################################################################")
    logger.debug(f"Running replication run from config file '{sim_config_file}'")
    with set_package_log_level(logging.INFO):
        simulators_and_read_filenames = run_simulator_from_sampler_per_channel_parallel(
            mk_run_dir=reads_dir, sim_params=sim_params, read_durations_per_channel=read_durations_per_channel, cycle_read_durations=(gap_sampler_type != "replication"),
            seq_end_time=sim_config.get("run_duration", None), use_nanosim_id=True, parallel_args={"n_jobs": -1},
        )
        # simulators_and_read_filenames = run_simulator_from_sampler_per_channel(
        #     mk_run_dir=reads_dir, sim_params=sim_params, read_durations_per_channel=read_durations_per_channel, cycle_read_durations=(gap_sampler_type != "replication"),
        #     seq_end_time=sim_config.get("run_duration", None), use_nanosim_id=True
        # )
    logger.debug(f"#################################################################")
    logger.debug(f"#################################################################")

    logger.info("Saving simulator statistics")
    write_simulator_stats([simulator for (simulator, _) in simulators_and_read_filenames], output_dir=run_dir)
    logger.info(f"Writing sequencing summary file '{seqsum_filename}'")
    convert_simfasta_dir_to_seqsum(reads_dir, seqsummary_filename=seqsum_filename)
    logger.info("Wrote sequencing summary file")

seqsum_filename = run_dir / "sequencing_summary.txt"
run_simulator(seqsum_filename)
create_figures(
    seqsum_filename, run_dir, delete_existing_figure_dir=True,
    group_column="all", ref_genome_path=None,
)

logger.debug(f"Done with replicating non-selective sequencing run script")
