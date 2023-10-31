"""
Plot an existing sequencing summary file (that possibly comes with a PAF file).

This is similar to the entrypoint `plot_seqsum`, but adapted to the cluster with sensible defaults.
It cleans up directories and plots figures from the simulator run as well, if they are present.
"""

import argparse
import logging
from pathlib import Path
from simreaduntil.seqsum_tools.seqsum_plotting import create_plots_for_seqsum
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, setup_logger_simple
from simreaduntil.shared_utils.plotting import ignore_tight_layout_warning
from simreaduntil.shared_utils.utils import delete_dir_if_exists, num_lines_in_file
from simreaduntil.simulator.utils import set_package_log_level
from simreaduntil.usecase_helpers.utils import create_figures, get_cleaned_seqsum_filename, remove_mux_scans_and_clean_if_inexistent

logger = setup_logger_simple(__name__)

add_comprehensive_stream_handler_to_logger(None)
logging.getLogger(__name__).setLevel(logging.DEBUG)

sequencing_summary_file = "sequencing_summary.txt"

paf_file = None
# paf_file = "~/ont_project_all/ont_project/runs/enrich_usecase/data/20190809_zymo_uncalled.paf"
# paf_file = "~/ont_project_all/ont_project/runs/enrich_usecase/data/deduped_20190809_zymo_uncalled.paf"

# sequencing_summary_file = Path("~/ont_project_all/ont_project/runs/data/20190809_zymo_seqsum.txt").expanduser()
# sim_run_dir = Path("zymo_realrun_plots/")

# import matplotlib.style as mplstyle
# import matplotlib
# matplotlib.use('pdf')
# mplstyle.use(["fast"])

with set_package_log_level(logging.DEBUG):
    seqsum_df = remove_mux_scans_and_clean_if_inexistent(sequencing_summary_file)
    
    create_figures(seqsum_df, run_dir=".", delete_existing_figure_dir=True, paf_file=paf_file, group_column="all" if paf_file is None else None)
    