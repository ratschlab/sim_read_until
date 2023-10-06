
from textwrap import dedent
from unittest.mock import MagicMock
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
from simreaduntil.shared_utils.dna import get_random_DNA_seq

from simreaduntil.usecase_helpers.utils import perfect_reads_gen, get_fasta_open_method, get_gen_nanosim_reads_cmd, normalize_fasta, random_nanosim_reads_gen, random_reads_gen

def test_read_gen(tmp_path):
    import itertools
    
    read_gen = random_nanosim_reads_gen(random_state=np.random.default_rng(2), length_range=(8, 17))
    list(itertools.islice(read_gen, 0, 5))
    
    ref_genome_path = tmp_path / "test.fa"
    with open(ref_genome_path, "w") as f:
        SeqIO.write([
            SeqIO.SeqRecord(seq=Seq(get_random_DNA_seq(100)), id="seq1"),
            SeqIO.SeqRecord(seq=Seq(get_random_DNA_seq(140)), id="seq2"),
        ], f, "fasta")
    read_gen = perfect_reads_gen(ref_genome_path, read_lens_range=(20, 30), random_state=np.random.default_rng(1))
    list(itertools.islice(read_gen, 0, 5))

def test_random_reads_gen():
    # just to check the function returns same thing when called twice with same random state, used for other tests
    import itertools
    get_random_reads = lambda: list(itertools.islice(random_reads_gen(random_state=np.random.default_rng(3)), 10))
    assert get_random_reads() == get_random_reads()

def test_gen_command():
    for use_slurm in [True, False]:
        get_gen_nanosim_reads_cmd("nanosim_dir", "nanosim_model_prefix", "ref_genome_path", "reads_dir", 100, n_procs=4, perfect=False, use_slurm=use_slurm)
        
def test_normalize_fasta(tmp_path, mocker):
    # write fake fasta file
    fasta_file = tmp_path / "test_in.fa"
    fasta_file_gz = tmp_path / "test_in.fa.gz"
    fasta_content = dedent(
        """\
        >seq1
        AGNcT
        >seq2
        AANNagAGCT
    """)
    
    with get_fasta_open_method(fasta_file, write=True)() as f:
        f.write(fasta_content)
    with get_fasta_open_method(fasta_file_gz, write=True)() as f:
        f.write(fasta_content)
    
    with get_fasta_open_method(fasta_file)() as f:
        assert str(f.read()) == fasta_content
    with get_fasta_open_method(fasta_file_gz)() as f:
        # f.read() not implemented
        assert "".join(line for line in f) == fasta_content
        
    out_fasta_file = tmp_path / "test_out.fa"
    out_fasta_file_gz = tmp_path / "test_out.fa.gz"
    
    mocked_state = mocker.patch("numpy.random.default_rng")
    mocked_state.choice.side_effect = ["C", "A", "G", "C", "A", "G", "A", "G", "A", "G"]
    normalize_fasta(fasta_file, out_fasta_file, random_state=mocked_state)
    normalize_fasta(fasta_file_gz, out_fasta_file_gz, random_state=mocked_state)
    
    normalized_fasta_content = dedent(
        """\
        >seq1 normalized
        AGCCT
        >seq2 normalized
        AAAGAGAGCT
    """)
    with get_fasta_open_method(out_fasta_file)() as f:
        assert str(f.read()) == normalized_fasta_content
    with get_fasta_open_method(out_fasta_file_gz)() as f:
        # f.read() not implemented
        assert "".join(line for line in f) == normalized_fasta_content
        
    # just one sequence
    normalize_fasta(fasta_file, out_fasta_file, seq_names=["seq2"], random_state=mocked_state)
    normalize_fasta(fasta_file_gz, out_fasta_file_gz, seq_names=["seq2"], random_state=mocked_state)
    normalized_fasta_content = dedent(
        """\
        >seq2 normalized
        AAAGAGAGCT
    """)
    with get_fasta_open_method(out_fasta_file)() as f:
        assert str(f.read()) == normalized_fasta_content
    with get_fasta_open_method(out_fasta_file_gz)() as f:
        # f.read() not implemented
        assert "".join(line for line in f) == normalized_fasta_content