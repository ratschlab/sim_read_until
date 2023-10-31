# superseed by seqsum_plotting.py
# todo: remove

import argparse
import logging
from pathlib import Path
import shutil

import toml
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, print_logging_levels, setup_logger_simple
from simreaduntil.shared_utils.utils import print_args
from simreaduntil.usecase_helpers.utils import plot_sim_stats

logger = setup_logger_simple(__name__)
"""module logger"""

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Plot simulator statistics. Use 'plot_seqsum' to plot sequencing summary statistics")
    parser.add_argument("run_dir", type=Path, help="Path to the run directory of the simulator with files ")
    parser.add_argument("--figure_dir", type=Path, help="Path to the figure directory, defaults to <run_dir>/figures, will overwrite files in there", default=None)
    parser.add_argument("--verbosity", type=str, help="log level for plotting commands", default="info")
    
    args = parser.parse_args(args)
    print_args(args, logger=logger)
    return args

def main():
    log_level = logging.DEBUG
    logging.getLogger(__name__).setLevel(log_level) # log level of this script (running as __main__)
    logging.getLogger("simreaduntil").setLevel(log_level) # log level of this script (running as __main__)
    add_comprehensive_stream_handler_to_logger(None, level=log_level)
    # print_logging_levels()

    args = parse_args()
    
    run_dir = args.run_dir
    assert run_dir.exists(), f"run_dir '{run_dir}' does not exist"
    figure_dir = args.figure_dir
    if figure_dir is None:
        figure_dir = run_dir / "figures"
    verbosity = args.verbosity
    
    figure_dir.mkdir(parents=True, exist_ok=True)
    
    logging.getLogger("simreaduntil").setLevel({"debug": logging.DEBUG, "info": logging.INFO, "warning": logging.WARNING, "error": logging.ERROR, "critical": logging.CRITICAL}[verbosity.lower()])

    plot_sim_stats(run_dir, figure_dir)
    
    logger.info("Done plotting simulator statistics")

if __name__ == "__main__":
    main()