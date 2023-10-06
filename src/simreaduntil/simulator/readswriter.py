"""
Write finished reads, e.g., to a file or list
"""

#%%
import logging
import os
from pathlib import Path
import shutil
import sys
import tempfile
from textwrap import indent
import threading
from Bio import SeqIO
from Bio.Seq import Seq
from simreaduntil.shared_utils.utils import is_empty_dir, setup_logger_simple

logger = setup_logger_simple(__name__)
"""module logger"""

#%%
class ReadsWriter:
    """
    Handle the writing of reads
    """
    def __init__(self):
        self._lock = threading.Lock()
        
    def write_read(self, read: SeqIO.SeqRecord):
        """
        write a read (or save it to write it later when flush is called)
        """
        with self._lock:
            self._write_read(read)
    
    def _write_read(self, read: SeqIO.SeqRecord):
        """thread-safe helper method for write_read"""
        raise NotImplementedError
    
    def flush(self):
        """Flush reads that were not yet written"""
        with self._lock:
            self._flush()
    
    def _flush(self):
        """thread-safe helper method for flush"""
        raise NotImplementedError

class SingleFileReadsWriter(ReadsWriter):
    """
    Write reads to one file (by default stdout), appending reads to the file as they are written (possibly with buffering)
    
    When pickling the file, the file is flushed. When reloading it, the filehandler is not restored and must be set directly.
    This class is useful for debugging by writing to sys.stdout.
    
    Args:
        reads_out_fh: filehandler to write to
        prefix: prefix to append to each id (useful if outputting to stdout)
    """
    def __init__(self, reads_out_fh=sys.stdout, prefix=None):
        super().__init__()
        
        self.fh = reads_out_fh
        self.prefix = '' if prefix is None else prefix
        
    # Flush reads, e.g. write outstanding reads to file by flushing the file handler
    def _flush(self):
        self.fh.flush()
        
    def __repr__(self):
        return f"SingleFileReadsWriter(filename='{self.fh.name}', prefix='{self.prefix}')"
        
    # write read to file
    def _write_read(self, read: SeqIO.SeqRecord):
        read.id = f"{self.prefix}{read.id}"
        SeqIO.write(read, self.fh, "fasta")
        
    # pickling support
    def __getstate__(self):
        if not self.fh.closed:
            # do not flush file if already closed
            self.flush()
        state = self.__dict__.copy()
        # pickling file handlers does not work as expected, e.g. when mode="w", reloading them will discard the existing file (which happens for example when running in parallel and serializing/deserializing the object)
        state["fh"] = None
        return state
    
class RotatingFileReadsWriter(ReadsWriter):
    """
    Write reads to a file, creating a new file whenever a maximum of reads is reached.
    
    First write reads to a temporary directory. Whenever a certain number of
    reads have been written, move the file to the output directory and start
    a new file. 
    It is essential to call .flush() before this object goes out of scope, otherwise 
    the reads will not be moved to the output directory. The context manager takes care of this.
    
    When pickling this object, the current file is flushed and moved to the output directory. Upon reloading,
    a new file will be opened.
    
    Args:
        output_dir: where to put completed read files; must not exist (since 
            files starting from zero idx will be written there)
            
        reads_out_prefix: prefix for filenames to write to
        max_reads_per_file: maximum reads to write per file
    """
    def __init__(self, output_dir, reads_out_prefix, max_reads_per_file=20):
        super().__init__()
        
        self.output_dir = Path(output_dir).absolute()
        # in case the directory is created by several channels, each writing to different files
        # check in this case that it is all empty
        if self.output_dir.exists():
            assert is_empty_dir(self.output_dir), f"reads output_dir '{self.output_dir}' must be empty"
        else:
            self.output_dir.mkdir()
        self.reads_out_prefix = reads_out_prefix
        self.max_reads_per_file = max_reads_per_file
        
        self.temp_dir = Path(tempfile.mkdtemp("readswriter")) # where to write files before moving them to output_dir
        self.current_file_idx = 0 # index of current file
        self.current_fh = None # current file handler
        self.current_nb_reads = 0 # number of reads written to current file
        
    def __repr__(self) -> str:
        return f"RotatingFileReadsWriter(output_dir='{self.output_dir}', reads_out_prefix='{self.reads_out_prefix}', max_reads_per_file={self.max_reads_per_file}, temp_dir='{self.temp_dir}', current_file_idx={self.current_file_idx}, current_nb_reads={self.current_nb_reads})"
    
    def _current_filename(self):
        return f"{self.reads_out_prefix}{self.current_file_idx}.fasta"
        
    def _open_new_file(self):
        assert self.current_fh is None
        filename = self.temp_dir / self._current_filename()
        logger.debug(f"Creating new reads file at location '{filename}'")
        self.current_fh = open(filename, "w")
        self.current_nb_reads = 0
    
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.flush()
        
    # close and move file to output_dir
    def _close_and_move_cur_file(self):
        if self.current_fh is None:
            return
        
        self.current_fh.close() # also flushes
        self.current_fh = None
        name = self._current_filename()
        assert not (self.output_dir / name).exists(), f"tried moving file from '{self.temp_dir / name}' to '{self.output_dir / name}', but destination already exists"
        shutil.move(self.temp_dir / name, self.output_dir / name) # os.rename does not work across file systems
        logger.debug(f"Moved file from '{self.temp_dir / name}' to '{self.output_dir / name}'")
        self.current_file_idx += 1
        
    def _flush(self):
        """
        Close current file and move it to output_dir
        
        The current file is closed and a new file is started because other tools reading the 
        FASTA only look for new files and don't check if existing files have been updated.
        
        It is essential to call this function before exiting, otherwise the reads will not be moved to the output directory.
        """
        self._close_and_move_cur_file()
        
    def _write_read(self, read: SeqIO.SeqRecord):
        """Write read to file"""
        
        # check if no current file -> open new file
        if self.current_fh is None:
            self._open_new_file()
        
        SeqIO.write(read, self.current_fh, "fasta")
        self.current_nb_reads += 1
        
        # check if enough reads written -> close current file
        if self.current_nb_reads >= self.max_reads_per_file:
            assert self.current_fh is not None
            self._close_and_move_cur_file()
    
    # pickling support
    def __getstate__(self):
        # pickling file handlers does not work as expected, e.g. when mode="w", reloading them will discard the existing file (which happens for example when running in parallel and serializing/deserializing the object)
        self.flush() # otherwise current_nb_reads not consistent, deletes current_fh
        assert self.current_fh is None
        state = self.__dict__.copy()
        return state
       
class ArrayReadsWriter(ReadsWriter):
    """
    Write reads to array, useful for debugging and testing
    """
    def __init__(self):
        super().__init__()
        
        self.reads = []
        self.output_dir = None # for compatibility with the ONTSimulator
        
    def __repr__(self) -> str:
        return f"ArrayReadsWriter(nb_reads={len(self.reads)})"
    
    def extended_repr(self):
        reads_str = "\n".join(f"{read_id} ({read_desc}): {read_seq}" for (read_id, read_seq, read_desc) in self.reads)
        return f"ArrayReadsWriter with reads\n{indent(reads_str, prefix='  ')}"
    
    def _flush(self):
        pass
    
    def _write_read(self, read: SeqIO.SeqRecord):
        self.reads.append((read.id, str(read.seq), read.description))

