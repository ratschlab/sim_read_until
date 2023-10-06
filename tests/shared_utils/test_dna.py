from textwrap import dedent
import pytest

from simreaduntil.shared_utils.dna import get_random_DNA_seq, get_nb_fasta_seqs, get_ref_lengths, get_reverse_complement

def test_get_random_DNA_seq():
    for l in [0, 2, 4, 7, 10]:
        assert len(get_random_DNA_seq(l)) == l

def test_get_reverse_complement():
    assert get_reverse_complement("ACGT") == "ACGT"
    assert get_reverse_complement("ACCGCAA") == "TTGCGGT"
    
def test_get_ref_lengths(tmp_path):
    fasta_example_path = tmp_path / "example.fasta"
    with open(fasta_example_path, "w") as f:
        print(dedent(f"""\
            >chr1 extra1
            AACTG
            >gen2-chr2 extra2
            ACTGAAA
        """), file=f)
    assert get_ref_lengths(fasta_example_path) == {"chr1": 5, "gen2-chr2": 7}
    assert get_nb_fasta_seqs(fasta_example_path) == 2