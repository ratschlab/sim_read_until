"""
Combines SimReadUntil with ReadFish

This script shows how to combine the SimReadUntil with ReadFish.

It first learns a gap sampler from an existing run and saves it.
Then, it runs the simulator in combination with ReadFish. on perfect reads generated from a reference genome, or from a reads file (if the config is adapted, e.g. NanoSim reads).
It creates the run results in the current directory, so switch to an appropriate directory before and use relative paths carefully in the config.

You can terminate the script early with `Ctrl+C` (Interrupt for the Jupyter kernel).
The script will print a warning if it cannot keep up with the acceleration factor. In this case, try to decrease the acceleration factor. The approach is real-time whenever the acceleration factor is at least 1.
You can find the reads in the appropriate directory.

You can try changing the acceleration factor, the frequency of the mux scans, the gap parameters, see the config file.
- perfect reads or the NanoSim reads generated above (which can also be perfect reads)
- write all reads to one file or in a rotating fashion
- load parameters extracted from an existing run (see below)

Run it with:
```
cd runs/enrich_usecase
source ~/ont_project_all/ont_project_venv/bin/activate
export PATH=~/ont_project_all/tools/bin:$PATH && which minimap2

mkdir dummy_random && cd dummy_random
ln -s ../data .
ln -s ../configs/dummy/random configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

cd ..
mkdir dummy_constant_gaps && cd dummy_constant_gaps
ln -s ../data .
ln -s ../configs/dummy/constant_gaps configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

cd ..
mkdir dummy_replication && cd dummy_replication
ln -s ../data .
ln -s ../configs/dummy/replication configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

cd ..
mkdir dummy_sampler_per_window && cd dummy_sampler_per_window
ln -s ../data .
ln -s ../configs/dummy/sampler_per_window configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py
# Press Ctrl+C to stop script, then plots are created

cd ..
mkdir full_run_sampler_per_window && cd full_run_sampler_per_window
ln -s ../data .
ln -s ../configs/full_run/sampler_per_window configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

cd ..
mkdir full_genome_run_sampler_per_window && cd full_genome_run_sampler_per_window
ln -s ../data .
ln -s ../configs/full_genome_run/sampler_per_window configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py

cd ..
mkdir chr202122_run && cd chr202122_run
ln -s ../data .
ln -s ../configs/chr202122_run/sampler_per_window configs
python ~/ont_project_all/ont_project/usecases/enrich_usecase.py
```
"""


import logging
from pathlib import Path
import sys
import warnings
import numpy as np
import toml

from simreaduntil.shared_utils.debugging_helpers import is_test_mode
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, print_logging_levels, setup_logger_simple
from simreaduntil.shared_utils.plotting import filter_seaborn_warnings
from simreaduntil.shared_utils.tee_stdouterr import TeeStdouterr
from simreaduntil.shared_utils.utils import delete_dir_if_exists, dill_dump, dill_load, print_cmd_and_run
from simreaduntil.simulator.utils import set_package_log_level
from simreaduntil.usecase_helpers import simulator_with_readfish
from simreaduntil.usecase_helpers.utils import create_simparams_if_inexistent, get_gap_sampler_method, plot_condor_log_file_metrics
from simreaduntil.usecase_helpers.utils import create_figures

logger = setup_logger_simple(__name__)

add_comprehensive_stream_handler_to_logger(None)
set_package_log_level(logging.INFO).__enter__()
print_logging_levels()
logging.getLogger(__name__).setLevel(logging.DEBUG)
# logging.getLogger().setLevel(logging.DEBUG) # warnings from everywhere, not desired

# import warnings
# warnings.filterwarnings("error")
filter_seaborn_warnings()

################################
## PARAMS
################################

if is_test_mode():
    # todo
    # import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/dummy_replication").expanduser())
    # import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/full_run_sampler_per_window/").expanduser())
    # import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/full_genome_run_sampler_per_window").expanduser())
    # import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/chr202122_run").expanduser())
    import os; os.chdir(Path("~/ont_project_all/ont_project/runs/enrich_usecase/full_genome_run_sampler_per_window").expanduser())
    # import os; os.chdir(Path("/Volumes/mmordig/ont_project/runs/enrich_usecase/full_genome_run_sampler_per_window").expanduser()) # not found due to samba bug not syncing config folder
sim_config_file = Path("configs/config.toml")
# ask_dir_deletion = True
ask_dir_deletion = False

################################
## SCRIPT
################################

sim_config = toml.load(sim_config_file)
run_dir = Path(sim_config["run_dir"])
# TeeStdouterr(run_dir / "stdouterr.txt").redirect()
logger.debug(f"Read in simulation config file '{sim_config_file}'")

ref_genome_path = sim_config.get("ref_genome_path", None)
sim_params_filename = Path(sim_config["sim_params_file"]) if "sim_params_file" in sim_config else None
seqsum_param_extr_file = Path(sim_config["seqsum_param_extr_file"]) if "seqsum_param_extr_file" in sim_config else None
gap_sampler_type = sim_config.get("gap_sampler_type", None)
gap_sampler_method = get_gap_sampler_method(gap_sampler_type, n_channels_full=sim_config["n_channels_full"])

logger.info("#"*80)
logger.info(f"Loaded config file '{sim_config_file}' with content:\n{sim_config_file.read_text()}")
logger.info("#"*80)
logger.info("#"*80)
logger.info(f"""Loading ReadFish config file with content:\n{Path(sim_config["readfish_config_file"]).read_text()}""")
logger.info("#"*80)

def create_minimap_index_if_inexistent():
    if sim_config["readfish_method"] != "unblock_all":
        readfish_config = toml.load(sim_config["readfish_config_file"])
        mmi_filename = readfish_config["conditions"]["reference"]
        if mmi_filename == "fake_mapper":
            logger.info(f"Skipping minimap2 index creation, using fake wrapper")
            return
        
        mmi_filename = Path(mmi_filename)
        if mmi_filename.exists():
            logger.info(f"Minimap2 index '{mmi_filename}' already exists, skipping minimap2 index creation")
        else:
            logger.debug(f"Creating minimap2 index at location '{mmi_filename}' for ReadFish from reference genome '{ref_genome_path}'")
            assert ref_genome_path is not None
            print_cmd_and_run(f"""minimap2 -d {mmi_filename} {ref_genome_path}""")
    else:
        logger.debug("Skipping minimap2 index (not needed)")

def run_readfish_simulation():
    logger.debug(f"#################################################################")
    logger.debug(f"#################################################################")
    logger.debug(f"Running the simulation from config file '{sim_config_file}' with ReadFish config file '{sim_config['readfish_config_file']}'")
    delete_dir_if_exists(run_dir, ask=ask_dir_deletion)
    # with set_package_log_level(logging.INFO):
    #     print_logging_levels()
    seqsum_file = simulator_with_readfish.main(sim_config_file)
    assert Path(seqsum_file).exists()
    logger.debug(f"#################################################################")
    logger.debug(f"#################################################################")
    
    return seqsum_file

# comment out as needed
create_minimap_index_if_inexistent() # comment this out if you want to use minimap2 to align to a reference
if sim_params_filename is None:
    logger.info(f"No sequencing summary specified to extract sim params, skipping sim params creation")
else:
    create_simparams_if_inexistent(sim_params_filename, seqsum_param_extr_file, sim_config.get("n_channels", None), gap_sampler_method)
assert run_dir / "sequencing_summary.txt" == run_readfish_simulation()

seqsum_filename = run_dir / "sequencing_summary.txt"
# read targets to enrich from file
figure_dir = run_dir / "figures"
delete_dir_if_exists(figure_dir, ask=ask_dir_deletion)
figure_dir.mkdir(exist_ok=True)

plot_condor_log_file_metrics(figure_dir)
create_figures(
    seqsum_filename, run_dir=run_dir, figure_dir=figure_dir,
    ref_genome_path=ref_genome_path, cov_thresholds=[1, 2, 3, 4],
    group_to_units={"target": toml.load(sim_config["readfish_config_file"])["conditions"]["0"]["targets"]},
)

logger.debug(f"Done with usecase script")