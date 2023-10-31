
import argparse
import itertools
import logging
import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
from simreaduntil.shared_utils.debugging_helpers import is_test_mode
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, print_logging_levels, setup_logger_simple
from simreaduntil.shared_utils.utils import print_args
from simreaduntil.simulator.readswriter import SingleFileReadsWriter
from simreaduntil.usecase_helpers.utils import random_nanosim_reads_gen

logger = setup_logger_simple(__name__)
"""module logger"""


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Generate dummy reads, only for testing purposes")
    parser.add_argument("reads_file", type=Path, help="Path to write reads to")
    parser.add_argument("--num_reads", type=int, help="Number of reads to generate", default=10_000)
    parser.add_argument("--length_range", type=str, help="Length range of reads", default="5_000,10_000")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing reads file")
    
    args = parser.parse_args(args)
    print_args(args, logger=logger)
    
    return args

def main():
    log_level = logging.DEBUG
    logging.getLogger(__name__).setLevel(log_level)
    add_comprehensive_stream_handler_to_logger(None, level=log_level)
    # print_logging_levels()
    
    if is_test_mode():
        os.chdir("server_client_cli_example")
        args = parse_args(args=["test.fasta", "--overwrite"]) #todo
    else:    
        args = parse_args()
    
    reads_file = args.reads_file
    overwrite = args.overwrite
    num_reads = args.num_reads
    assert num_reads > 0, f"num_reads {num_reads} must be > 0"
    length_range = args.length_range
    length_range = tuple(map(int, length_range.split(",")))
    assert len(length_range) == 2, f"length_range {length_range} must have length 2"
    assert length_range[0] < length_range[1], f"length_range {length_range} must be increasing"
    
    assert reads_file.suffix == ".fasta", f"reads_file '{reads_file}' must end with '.fasta'"
    if reads_file.exists():
        if overwrite:
            # logger.warning(f"Overwriting existing reads file '{reads_file}'")
            reads_file.unlink()
        else:
            raise FileExistsError(f"reads_file '{reads_file}' already exists, use --overwrite to overwrite")
    
    reads_gen = random_nanosim_reads_gen(length_range=length_range)
    with open(reads_file, "w") as fh:
        writer = SingleFileReadsWriter(fh)
        for (read_id, seq) in tqdm(itertools.islice(reads_gen, num_reads), total=num_reads, desc="Writing read: "):
            writer.write_read(SeqIO.SeqRecord(id=str(read_id), seq=Seq(seq)))
    logger.info(f"Done writing reads to file '{reads_file}'")

if __name__ == "__main__":
    main()