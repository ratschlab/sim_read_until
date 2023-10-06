import pytest
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from simreaduntil.shared_utils.dna import get_random_DNA_seq
from simreaduntil.simulator.readpool import ReadPoolFromIterable, ReadPoolFromIterablePerChannel, ReadPoolFromFile, reads_from_file_gen, NoReadLeft

gen_from_list = pytest.helpers.gen_from_list


@pytest.fixture
def dummy_reads_fasta(tmp_path):
    """
    Creates a dummy fasta file with reads
    """
    filename = tmp_path / "test111.fasta"
    with open(filename, "w") as f:
        SeqIO.write(SeqIO.SeqRecord(Seq("AAACCTGG"), id="read1"), f, "fasta")
        SeqIO.write(SeqIO.SeqRecord(Seq("CCCCCGGTT"), id="read2"), f, "fasta")
        return filename
        
        
def test_ReadPoolFromIterable():
    # read sequences from list/generator
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "GGGAAATCCGAA"), ("read2", "AAACCTGGTTAGG"))))
    assert read_pool.nb_reads_returned == 0
    assert read_pool.get_new_read() == ("read1", "GGGAAATCCGAA")
    assert read_pool.nb_reads_returned == 1
    assert read_pool.get_new_read() == ("read2", "AAACCTGGTTAGG")
    assert read_pool.nb_reads_returned == 2
    with pytest.raises(NoReadLeft):
        read_pool.get_new_read()
    assert read_pool.nb_reads_returned == 2
        
def test_ReadPoolFromIterablePerChannel():
    random_state = np.random.default_rng(2)
    random_seqs_gen_from_lens = lambda lengths: (get_random_DNA_seq(n=n, random_state=random_state) for n in lengths)
    
    read_pool = ReadPoolFromIterablePerChannel({1: random_seqs_gen_from_lens([2, 3, 5]), 2: random_seqs_gen_from_lens([8, 9])})
    with pytest.raises(KeyError):
        read_pool.get_new_read()
        
    assert len(read_pool.get_new_read(1)) == 2
    assert len(read_pool.get_new_read(2)) == 8
    assert len(read_pool.get_new_read(1)) == 3
    assert len(read_pool.get_new_read(1)) == 5
    with pytest.raises(NoReadLeft):
        read_pool.get_new_read(1)
    # do twice
    with pytest.raises(NoReadLeft):
        read_pool.get_new_read(1)
        
    assert len(read_pool.get_new_read(2)) == 9
    with pytest.raises(NoReadLeft):
        read_pool.get_new_read(2)
    
def test_reads_from_file_gen(dummy_reads_fasta):
    # no shuffle
    assert list(reads_from_file_gen(dummy_reads_fasta)) == [('read1', 'AAACCTGG'), ('read2', 'CCCCCGGTT')]
    
    random_state = np.random.default_rng(2)
    assert list(reads_from_file_gen(dummy_reads_fasta, shuffle_rand_state=random_state)) == [('read1', 'AAACCTGG'), ('read2', 'CCCCCGGTT')]
    random_state = np.random.default_rng(5)
    assert list(reads_from_file_gen(dummy_reads_fasta, shuffle_rand_state=random_state)) == [('read2', 'CCCCCGGTT'), ('read1', 'AAACCTGG')]
    
def test_ReadPoolFromFile(dummy_reads_fasta):
    read_pool = ReadPoolFromFile(dummy_reads_fasta)
    assert not read_pool.definitely_empty
    
    print(read_pool)
    
    assert read_pool.get_new_read() == ("read1", "AAACCTGG")
    assert read_pool.get_new_read() == ("read2", "CCCCCGGTT")
    with pytest.raises(NoReadLeft):
        read_pool.get_new_read()
    assert read_pool.definitely_empty
    
