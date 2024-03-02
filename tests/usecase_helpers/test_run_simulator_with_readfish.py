

import logging
import os
from matplotlib import pyplot as plt
import pandas as pd
import toml
from pathlib import Path
from simreaduntil.shared_utils.utils import print_cmd_and_run
from simreaduntil.simulator.simulator import plot_sim_actions
from simreaduntil.simulator.utils import set_package_log_level

from simreaduntil.usecase_helpers.simulator_with_readfish import main as run_readfish

def create_minimap_index_if_inexistent(mmi_filename, ref_genome_path):
    mmi_filename = Path(mmi_filename)
    if mmi_filename.exists():
        return False
    else:
        print_cmd_and_run(f"""minimap2 -d {mmi_filename} {ref_genome_path}""", verbose=False)
        return True
    
def test_simulator_with_readfish(shared_datadir, tmp_path):
    # we run the simulation for 10 seconds
        
    # tmp_path already has a data dir due to shared_datadir, so we use a subdirectory
    run_dir = tmp_path / "run_dir"
    run_dir.symlink_to(shared_datadir / "run_dir", target_is_directory=True)
    os.chdir(run_dir)
    
    create_minimap_index_if_inexistent("data/chm13v2.0_normalized1000000firsttwo.mmi", ref_genome_path="data/chm13v2.0_normalized1000000firsttwo.fa.gz")
    with set_package_log_level(logging.WARNING):
        assert run_readfish(Path("configs/config.toml")).resolve() == Path("simulator_run/sequencing_summary.txt").resolve()
    
    assert Path("simulator_run/reads").exists()
    assert Path("simulator_run/sequencing_summary.txt").exists()
    assert Path("simulator_run/live_sequencing_summary.txt").exists()
    
    assert Path("simulator_run/sequencing_summary.txt").read_text() == Path("simulator_run/live_sequencing_summary.txt").read_text()
    
    action_results_df = pd.read_csv("simulator_run/action_results.csv", sep="\t")
    plot_sim_actions(action_results_df, close_figures=True)
    