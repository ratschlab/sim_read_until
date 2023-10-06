"""
Common helper functions for dealing with DNA
"""

from typing import Dict
import numpy as np
import pysam

from simreaduntil.shared_utils.utils import get_file_content, is_empty_file

def get_random_DNA_seq(n, random_state=np.random.default_rng(2)):
    """
    Random DNA sequence of length n
    
    Args:
        n: length of the sequence
        random_state: random state; if not provided, repeated calls to this function return different values
    """
    return "".join(random_state.choice(list("ACGT"), n, replace=True))

def get_reverse_complement(seq):
    """Reverse complement of DNA sequence"""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def get_ref_lengths(ref_genome_path) -> Dict[str, int]:
    """
    Get chromosome lengths in a dict
    
    It supports .gz file
    """
    if is_empty_file(ref_genome_path):
        return {}
    with pysam.FastaFile(ref_genome_path) as fasta:
        return {ref: fasta.get_reference_length(ref) for ref in fasta.references}
    
def get_nb_fasta_seqs(fasta_filename):
    """
    Get number of sequences in a fasta file
    
    It supports .gz file
    """
    if is_empty_file(fasta_filename):
        return 0
    
    try:
        with pysam.FastaFile(fasta_filename) as fasta:
            return fasta.nreferences
    except OSError:
        # for debugging
        # pysam raises OSError for empty files
        print(f"File content of '{fasta_filename}':\n{get_file_content(fasta_filename)}EOF")
        raise

def get_sequence_names_in_fasta(fasta_filename):
    """
    Get sequence names in a fasta file

    Args:
        fasta_filename: filename of fasta file

    Returns:
        sequence names in FASTA
    """
    if is_empty_file(fasta_filename):
        return []
    
    # first time (30s), then 0s (due to cache)
    with pysam.FastaFile(fasta_filename) as fasta_file:
        return fasta_file.references

    # pyfastx similar to pyfasta, slightly larger index, but full sequence names
    # import pyfastx
    # for seq in pyfastx.Fasta(fasta_filename):
    #     seq.name

    # plain file reading (58s)
    # sequence_names = []
    # with gzip.open(fasta_filename, "rt") as fasta_file:
    #     for line in tqdm.tqdm(fasta_file):
    #         if line.startswith(">"):
    #             sequence_names.append(line[1:].strip())
    # return sequence_names

    # using BioSeq is slower (1m13s), but gives access to the full sequence names (record.description)
    # from Bio import SeqIO
    # sequence_names = []
    # with get_fasta_open_method(fasta_filename)() as f:
    #     for record in tqdm.tqdm(SeqIO.parse(f, "fasta")):
    #         sequence_names.append(record.id)
    # return sequence_names
