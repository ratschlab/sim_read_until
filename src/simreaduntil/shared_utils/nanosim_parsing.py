"""
Helper functions for parsing NanoSim read ids, as well as some functions modified from NanoSim
"""

# from Bio import SeqIO
import numpy as np

class NanoSimId:
    """
    Represent a NanoSim id, for aligned and perfect reads
    
    Integer-like strings will be converted.
    """
    def __init__(self, chrom, ref_pos, read_nb, direction, ref_len, head_len=0, tail_len=0, read_type="aligned"):
        assert read_type in ["aligned", "perfect"]
        assert direction in ["R", "F"]
        ref_pos = int(ref_pos)
        ref_len = int(ref_len)
        
        head_len = int(head_len)
        tail_len = int(tail_len)
        # if read_type == "perfect":
        assert head_len == 0
        assert tail_len == 0
        
        self.chrom = chrom
        self.ref_pos = ref_pos
        self.read_nb = read_nb # may not be an integer
        self.direction = direction
        self.ref_len = ref_len
        self.head_len = head_len
        self.tail_len = tail_len
        self.read_type = read_type
        
    @staticmethod
    def from_str(read_id: str):
        """
        Parse NanoSim non-chimeric read id, for aligned and perfect reads
        
        Note: Perfect reads are not chimeric (in NanoSim).
        
        E.g.
        chr2_920875_perfect_proc0:1_F_0_8346_0
        chr2_649870_aligned_proc0:2_F_0_10399_0
        
        This is of the form {genome}-{chromosome}_{ref_pos}_{read_type}_{read_nb}_{strand}_{head_len}_{segment_lengths}_{tail_len}
        The ref_pos is 0-based and the read spans [ref_pos:ref_pos+ref_len] on the forward strand, independent of the the direction which is F or R (forward, reverse).
        We assume here that the head and tail flanking lengths are 0.
        read_nb: proc{process_nb}:{read_nb_for_process}
        
        Note (not relevant here): for chimeric reads, '{genome-chromosome}_{position}' and 'segment_lengths' are joined by ";" etc.
        """
        chrom, ref_pos, read_type, read_nb, direction, head_len, ref_len, tail_len = read_id.split("_")
        return NanoSimId(chrom, ref_pos, read_nb, direction, ref_len, head_len, tail_len, read_type)
        
    @staticmethod
    def is_valid(read_id: str):
        """
        Simple check whether something is a NanoSim id, has false positives
        """
        return len(read_id.split("_")) == 8
    
    def change_ref_len(self, new_ref_len: int):
        """
        Changes the reference length of the read, i.e. ref_pos and ref_len depending on strandness, and adds a "m" to the read_nb, in-place.
        
        The read is still uniquely identifiable because a unique read number is encoded in self.read_nb.
        
        Args:
            new_ref_len: new reference length, must be less or equal to the current reference length
                the method returns early if the reference length is already the same.
                
        Returns:
            self or None
        """
        assert self.ref_len >= new_ref_len >= 0, f"new length {new_ref_len} should be in range [0, {self.ref_len}]"
        
        if self.ref_len == new_ref_len:
            return self
        
        # modify the start location which is with respect to the forward strand
        if self.direction == "R":
            self.ref_pos += self.ref_len - new_ref_len
        self.ref_len = new_ref_len
        self.read_nb += "m" # to indicate that the read has been modified
        
        return self
    
    def __repr__(self):
        return f"{self.chrom}_{self.ref_pos}_{self.read_type}_{self.read_nb}_{self.direction}_{self.head_len}_{self.ref_len}_{self.tail_len}"


def normalize_seq_name(seq_name):
    """
    Normalize sequence name (or chromosome name) from FASTA record
    
    Adapted from Nanosim
    
    Args:
        seq_name: a fasta header
    Returns:
        Split by whitespaces and underscore and join the fragments by "-"; return everything before the first dot "."
    """
    import re
    info = re.split(r'[_\s]\s*', seq_name) # \s = whitespace character
    chr_name = "-".join(info)
    return chr_name.split(".")[0]

def case_convert_dna(seq, random_state=np.random.default_rng(2)):
    """
    Translate lower to upper and replace N by any base, Y, D also
    
    Adapted from Nanosim function "case_convert"
    
    Args:
        seq: DNA sequence
        random_state: random state; if not provided, repeated calls to this function return different values
        
    Returns:
        DNA sequence with upper case letters and N, Y, D replaced by a random base
    """
    base_code = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
                 'N': ['A', 'T', 'C', 'G'], 'X': ['A', 'T', 'C', 'G']}

    up_string = seq.upper()
    up_list = list(up_string)
    for i in range(len(up_list)):
        if up_list[i] in base_code:
            up_list[i] = random_state.choice(base_code[up_list[i]])
    out_seq = ''.join(up_list)

    return out_seq
