from pathlib import Path
from unittest.mock import MagicMock
import pytest
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from simreaduntil.shared_utils.dna import get_random_DNA_seq
from simreaduntil.simulator.readpool import ReadPoolFromIterable, ReadPoolFromIterablePerChannel, ReadPoolFromFile, ThreadedReadPoolWrapper, reads_from_file_gen, NoReadLeftException

gen_from_list = pytest.helpers.gen_from_list


@pytest.fixture
def dummy_reads_fastas(tmp_path):
    """
    Creates a dummy fasta file with reads
    """
    filename1 = tmp_path / "test1.fasta"
    with open(filename1, "w") as f:
        SeqIO.write(SeqIO.SeqRecord(Seq("AAACCTGG"), id="read1"), f, "fasta")
        SeqIO.write(SeqIO.SeqRecord(Seq("CCCCCGGTT"), id="read2"), f, "fasta")

    filename2 = tmp_path / "test2.fasta"
    with open(filename2, "w") as f:
        SeqIO.write(SeqIO.SeqRecord(Seq("AAACCTGG"), id="read3"), f, "fasta")
        SeqIO.write(SeqIO.SeqRecord(Seq("CCCCCGGTT"), id="read4"), f, "fasta")
    
    return filename1, filename2
        
        
def test_ReadPoolFromIterable():
    # read sequences from list/generator
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "GGGAAATCCGAA"), ("read2", "AAACCTGGTTAGG"))))
    assert read_pool.nb_reads_returned == 0
    assert read_pool.get_new_read() == ("read1", "GGGAAATCCGAA")
    assert read_pool.nb_reads_returned == 1
    assert read_pool.get_new_read() == ("read2", "AAACCTGGTTAGG")
    assert read_pool.nb_reads_returned == 2
    with pytest.raises(NoReadLeftException):
        read_pool.get_new_read()
    assert read_pool.nb_reads_returned == 2
    
    read_pool.finish()
        
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
    with pytest.raises(NoReadLeftException):
        read_pool.get_new_read(1)
    # do twice
    with pytest.raises(NoReadLeftException):
        read_pool.get_new_read(1)
        
    assert len(read_pool.get_new_read(2)) == 9
    with pytest.raises(NoReadLeftException):
        read_pool.get_new_read(2)
    
    read_pool.finish()
    
def test_reads_from_file_gen(dummy_reads_fastas):
    # no shuffle
    assert list(reads_from_file_gen(dummy_reads_fastas[0])) == [('read1', 'AAACCTGG'), ('read2', 'CCCCCGGTT')]
    
    random_state = np.random.default_rng(2)
    assert list(reads_from_file_gen(dummy_reads_fastas[0], shuffle_rand_state=random_state)) == [('read1', 'AAACCTGG'), ('read2', 'CCCCCGGTT')]
    random_state = np.random.default_rng(5)
    assert list(reads_from_file_gen(dummy_reads_fastas[0], shuffle_rand_state=random_state)) == [('read2', 'CCCCCGGTT'), ('read1', 'AAACCTGG')]
    
def test_ReadPoolFromFile(dummy_reads_fastas):
    read_pool = ReadPoolFromFile(dummy_reads_fastas[0])
    assert not read_pool.definitely_empty
    
    print(read_pool)
    
    assert read_pool.get_new_read() == ("read1", "AAACCTGG")
    assert read_pool.get_new_read() == ("read2", "CCCCCGGTT")
    with pytest.raises(NoReadLeftException):
        read_pool.get_new_read()
    assert read_pool.definitely_empty
    
def test_ReadPoolFromFileThreaded(dummy_reads_fastas):
    reads_dir = Path(dummy_reads_fastas[0]).parent
    
    for (file_obj, num_reads) in [(dummy_reads_fastas[0], 2), (reads_dir, 4)]:
    
        assert ReadPoolFromFile.can_handle(file_obj)
    
        read_pool = ThreadedReadPoolWrapper(ReadPoolFromFile(file_obj), queue_size=3)
        repr(read_pool)
        
        [read_pool.get_new_read() for _ in range(num_reads)]
        with pytest.raises(NoReadLeftException):
            read_pool.get_new_read()
        with pytest.raises(NoReadLeftException):
            read_pool.get_new_read()
            
        assert read_pool.definitely_empty
        
        read_pool.finish()
    
def test_ReadPoolFromFileThreaded_NonEmpty(dummy_reads_fastas):
    reads_dir = Path(dummy_reads_fastas[0]).parent
    read_pool = ThreadedReadPoolWrapper(ReadPoolFromFile(reads_dir), queue_size=2)
    read_pool.get_new_read()
    
    read_pool.finish() # joins the thread
    # there should still be one read left, check that the read pool thread terminates properly
    
def pyslow5_is_available():
    try:
        import pyslow5
        return True
    except ImportError:
        return False
    
# @pytest.mark.skipif(not pyslow5_is_available(), reason="pyslow5 is not installed")
# def test_Slow5ReadPool(mocker, tmp_path):
#     import tempfile
#     mock = MagicMock() # cannot use Mock() because it doesn't have __iter__
#     mocker.patch("simreaduntil.simulator.readpool.get_slow5_reads_gen", return_value=mock)
#     mock.__iter__.return_value = gen_from_list((("read1", [1, 2, 3]), ("read2", [4, 5])))
#     slow5_dir = tmp_path / "slow5_dummy"
#     slow5_dir.mkdir()
#     (slow5_dir / "test.slow5").touch() # so it iterates over one file
    
#     assert Slow5ReadPool.can_handle(slow5_dir)
    
#     read_pool = Slow5ReadPool(slow5_dir, 1)
#     assert read_pool.get_new_read() == ("read1", [1, 2, 3])
#     assert read_pool.get_new_read() == ("read2", [4, 5])
#     with pytest.raises(NoReadLeftException):
#         read_pool.get_new_read()
#     assert read_pool.definitely_empty