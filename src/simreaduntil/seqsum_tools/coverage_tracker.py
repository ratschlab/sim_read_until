"""
Track coverage over time by parsing the NanoSim read ids
"""

import functools
from numbers import Number
import os
from typing import Dict, Iterable, List, Optional, Union, Tuple, Any
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
import tqdm
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId, normalize_seq_name
from simreaduntil.shared_utils.plotting import get_fig_ax

logger = setup_logger_simple(__name__)
"""module logger"""

def get_fraction_at_least(arr, threshold, length=None):
    """
    Compute the fraction of values in an array that are at least a given threshold
    
    Args:
        arr: array of values
        threshold: threshold
        length: if given, the length of the array to divide by, else len(arr)
    
    Returns:
        fraction of values in arr that are at least threshold
    """
    if length is None:
        return (arr >= threshold).mean()
    return (arr >= threshold).sum() / length

def plot_cov_per_position_lineplot(cov_per_position, n_points=100, ax=None, target_coverage: Optional[int]=None, cov_mean=None, cov_median=None, blocksize=1, **plot_args) -> plt.Axes:
    """
    Plot coverage of each basepair, with at most n_points
    
    Positions start at 0.
    
    Args:
        cov_per_position: coverage of each basepair, in the order of genome position
        n_points: maximum number of points to plot
        ax: axes to plot on
        target_coverage: if not None, plot a horizontal line at this coverage as a red dotted line
        cov_mean: mean coverage to plot, if None, mean of cov_per_position
        cov_median: median coverage to plot, if None, median of cov_per_position
        blocksize: blocksize that was used to compute cov_per_position, one position is one block
        plot_args: passed to sns.lineplot
        
    Returns:
        axis
    """
    fig, ax = get_fig_ax(ax, **plot_args)
    
    if blocksize > 1:
        assert cov_mean is not None
        assert cov_median is not None
    
    every = max(1, round(len(cov_per_position)/n_points))
    x_indices = np.arange(0, len(cov_per_position), every)
    
    cov_per_position = np.array(cov_per_position)
    sns.lineplot(x=x_indices * blocksize, y=cov_per_position[x_indices], ax=ax, marker="o")
    cov_mean = cov_per_position.mean() if cov_mean is None else cov_mean
    ax.axhline(y=cov_mean, linestyle="dotted", color="black", label=f"mean cov {cov_mean:.2f}")
    cov_median = np.median(cov_per_position) if cov_median is None else cov_median
    ax.axhline(y=cov_median, linestyle="dashed", color="black", label=f"median cov {cov_median:.2f}")
    if target_coverage is not None:
        ax.axhline(y=target_coverage, linestyle="dotted", color="red", label="target cov")
    ax.set_title(f"Coverage of {'basepairs' if blocksize == 1 else 'blocks'}")
    ax.set_xlabel('Position' if blocksize == 1 else 'Block start')
    ax.set_ylabel('Coverage')
    ax.set_xlim(0)
    # ax.set_ylim(0)
    ax.legend()

    return ax

def _compute_bin_edges(vals, discrete):
    """
    Compute bin edges for histogram if discrete is True, else return None
    Similar to seaborn's 'discrete' argument in histplot
    """
    if discrete:
        # vals = vals.dropna()
        start, stop = vals.min(), vals.max()
        return np.arange(start - .5, stop + 1.5)
    return None

def plot_covs_hist(cov_per_position, ax=None, target_coverage: Optional[int]=None, cov_mean=None, cov_median=None, **plot_args):
    """
    Plot coverage distribution
    The mean is plotted in dotted black, the median in dashed black
    
    Args:
        cov_per_position: coverage array
        ax: matplotlib axes to plot on
        target_coverage: if not None, plot a vertical line at this coverage as a red line
        cov_mean: mean coverage to plot, if None, mean of cov_per_position
        cov_median: median coverage to plot, if None, median of cov_per_position
        plot_args: passed to get_fig_ax
        
    Returns:
        axis
    """
    fig, ax = get_fig_ax(ax, **plot_args)
    if len(cov_per_position) > 0:
        discrete = True if (max(cov_per_position) < 100) else False
        # seaborn much slower than plt, see https://github.com/mwaskom/seaborn/issues/3325 
        # sns.histplot(covs, discrete=discrete, ax=ax)
        # sns.histplot(covs, bins=range(0, 12, 1), ax=ax)
        # ax.hist(covs, bins=range(0, 12, 1))
        cov_per_position = np.array(cov_per_position)
        ax.hist(cov_per_position, bins=_compute_bin_edges(cov_per_position, discrete))
        cov_mean = cov_per_position.mean() if cov_mean is None else cov_mean
        ax.axvline(x=cov_mean, linestyle="dotted", color="black", label=f"mean cov {cov_mean:.2f}")
        cov_median = np.median(cov_per_position) if cov_median is None else cov_median
        ax.axvline(x=cov_median, linestyle="dashed", color="black", label=f"median cov {cov_median:.2f}")
    if target_coverage is not None:
        ax.axvline(x=target_coverage, linestyle="dotted", color="red", label="target cov")
    ax.set_title("Coverage distribution")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Number")
    ax.legend()

    return ax

class BasepairCoverageTracker:
    """
    Track coverage per chromosome in a genome or metagenome, given the read position (e.g. parsed from NanoSim id or from a PAF file)
    
    Make sure that the counts per bp do not exceed np.uint8, i.e. < 255.
    """
    def __init__(self, coverage_per_chrom):
        self.coverage_per_chrom = coverage_per_chrom
        
    @classmethod
    def empty_from_ref_genomes(cls, ref_genomes: Union[Dict[str, os.PathLike], os.PathLike], **kwargs):
        """
        Create empty coverage arrays (with shapes of the reference genomes) in the file
        
        For a dict {genome_name: genome_file}, it prefixes each chromosome name with the genome name.
        """
        coverage_per_chrom = {}

        if isinstance(ref_genomes, dict):
            # metagenome mode
            for genome_name, genome_path in ref_genomes.items():
                coverage_per_chrom.update(cls._zero_arr_per_chrom(genome_path, seq_prefix=f"{genome_name}-"))
        else:
            # genome mode: no seq prefix
            coverage_per_chrom = cls._zero_arr_per_chrom(ref_genomes)
        
        return cls(coverage_per_chrom=coverage_per_chrom, **kwargs)
    
    @classmethod
    def empty_from_lens(cls, chrom_lens: Dict[str, int], **kwargs):
        """
        Create empty coverage arrays according to lengths
        """
        coverage_per_chrom = {genome_name: np.zeros(genome_len, dtype=np.uint8) for (genome_name, genome_len) in chrom_lens.items()}
        return cls(coverage_per_chrom=coverage_per_chrom, **kwargs)
    
    @staticmethod
    def _zero_arr_per_chrom(genome_filename, seq_prefix="") -> Dict[str, np.ndarray]:
        """
        Initialize an array with zeros for each chromosome in the genome
        
        Args:
            genome_filename: fasta file to identify chromosome names
            seq_prefix: prefix to add to each chromosome name
            
        Returns:
            dict mapping chromosome to a zero array
        """
        with pysam.FastaFile(genome_filename) as fasta:
            # use uint8 because counts should be less than 256 usually
            return {seq_prefix + normalize_seq_name(id): np.zeros(fasta.get_reference_length(id), dtype=np.uint8) for id in fasta.references}

    def get_chrom_start_len(self, read_id) -> Optional[Tuple[Any, int, int]]:
        """
        Map a read to its chromosome, start and length on the reference
        
        Args:
            read_id: a NanoSim read id
        
        Returns:
            Tuple (chrom, ref_start, ref_len) with respect to forward strand; None if could not be mapped
        """
        raise NotImplementedError()
    
    def add_reads(self, read_ids: Iterable, verbose=False):
        """
        Add new reads that were sequenced (as this influences coverage) 
        
        Args:
            read_ids: list or generator of read ids, should be parsable by NanoSimId.parse, 
                ref_len should be adapted if read was rejected in the middle
                
            verbose: whether to output regularly how many reads have been added
        
        Returns: 
            per chromosome, (number of reads added, number of basepairs added)
        """
        nb_reads = {chrom: 0 for chrom in self.coverage_per_chrom.keys()}
        nb_basepairs = {chrom: 0 for chrom in self.coverage_per_chrom.keys()}
        for read_id in tqdm.tqdm(read_ids) if verbose else read_ids:
            mapping_res = self.add_read(read_id)
            if mapping_res is not None:
                chrom, ref_start, ref_len = mapping_res
                nb_reads[chrom] += 1
                nb_basepairs[chrom] += ref_len
        
        if verbose:
            logger.info(f"Added {sum(nb_reads.values())} reads with {sum(nb_basepairs.values())} bps")
            
        return nb_reads, nb_basepairs
    
    def add_read(self, read_id):
        """Add a single read with this read_id
        """
        mapping_res = self.get_chrom_start_len(read_id)
        if mapping_res is None:
            logger.debug(f"Could not map read with id '{read_id}'")
            return None
        chrom, ref_start, ref_len = mapping_res
        # note: assumes that ref_len is read reference length, i.e. if read was rejected in the middle, the ref length was adapted accordingly
        assert ref_start + ref_len <= len(self.coverage_per_chrom[chrom]), f"Read {read_id} mapped to {chrom}:{ref_start}-{ref_start+ref_len} but chromosome length is {len(self.coverage_per_chrom[chrom])}, circular genome currently not supported"
        self.coverage_per_chrom[chrom][ref_start:ref_start+ref_len] += 1
        return chrom, ref_start, ref_len

    def avg_chrom_coverage(self, chrom) -> float:
        """
        Compute average coverage of a chromosome
        
        Args:
            chrom: chromosome
            
        Returns:
            average coverage for the chromosome
        """
        return self.coverage_per_chrom[chrom].mean()
    
    def get_fraction_cov_atleast(self, threshold, chroms: Optional[List]=None) -> float:
        """
        Compute the fraction of basepairs above a given threshold
        
        Args:
            threshold: coverage threshold
            chroms: chromosomes to consider
        """
        if chroms is None:
            chroms = list(self.coverage_per_chrom.keys())
        if len(chroms) == 0:
            logger.warning("get_fraction_cov_atleast called with empty chroms, returning 1.0")
            return 1.0
        return sum((self.coverage_per_chrom[chrom] >= threshold).sum(dtype=np.uint64) for chrom in chroms) / sum(len(self.coverage_per_chrom[chrom]) for chrom in chroms)
        
    def get_chrom_lens(self) -> Dict[str, int]:
        """
        Get lengths of chromosomes
        """
        return {chrom: len(arr) for (chrom, arr) in self.coverage_per_chrom.items()}
    
    # mean coverage is shown in red
    # if target coverage is given, it is plotted as a dashed line
    def plot_state(self, plot_type, target_coverage=None, **kwargs):
        """
        Plot coverage per chromosome with one axis per chromosome, as a line plot (over chromosome position) or histogram
        
        Args:
            plot_type: "line" or "hist"
            target_coverage: if given, plot a dashed line at this coverage level
            kwargs: passed to plot_cov_per_position_lineplot or plot_covs_hist
        """
        assert plot_type in ["line", "hist"]
        
        import matplotlib.pyplot as plt
        
        num_chroms = len(self.coverage_per_chrom)
        fig, axes = plt.subplots(ncols=num_chroms, figsize=(6, 1 + 3*num_chroms), squeeze=False)
        for ((chrom, covs_per_bp), ax) in zip(self.coverage_per_chrom.items(), axes.flat):
            plot_fcn = plot_cov_per_position_lineplot if plot_type == "line" else plot_covs_hist
            plot_fcn(covs_per_bp, ax=ax, target_coverage=target_coverage, **kwargs)
            if target_coverage is None:
                ax.set_title(chrom)
            else:
                fraction_at_least = get_fraction_at_least(covs_per_bp, target_coverage)
                ax.set_title(f"{chrom}: {fraction_at_least*100:.2f}%")
            
        [fig.delaxes(ax) for ax in axes.flat[num_chroms:]] # remove extra axes
        fig.suptitle("Coverage per basepair")
        
        return fig
    
def _add_count_to_blocked_array(arr, start, end, blocksize):
    """
    Increment count of region [start, end) in a block array
    
    A block array keeps track of counts per block of size blocksize, the last block possibly being shorter.
    Each block is represented by a single integer and counts the number of basepairs covered by a read in that block.
    """
    assert 0 <= start <= end <= len(arr) * blocksize
    start_block = start // blocksize
    end_block = end // blocksize # last block
    arr[start_block:end_block] += blocksize
    # correct start block
    arr[start_block] -= start % blocksize
    if end_block < len(arr):
        # add end block
        arr[end_block] += end % blocksize

def _divide_round_up(a: int, b: int):
    # divide a by b, rounding up
    return (a - 1) // b + 1

class BlockwiseCoverageTracker:
    """
    Track coverage per chromosome in a genome or metagenome in blocks, given the read position (e.g. parsed from NanoSim id or from a PAF file)
    
    Each chromosome [0, L) is divided into blocks of size L (called blocksize), the last possibly being shorter.
    When we add a read that overlaps with m basepairs of a block, we increment that block's count by m.
    
    For larger blocksize, the metric get_fraction_cov_atleast becomes less representative of the basepairs covered at least x times since
    it requires the average of a block to be above a threshold to count all basepairs in that block.
    Larger blocksize reduces memory usage and increases speed.
    
    Make sure that the number of basepairs covered by reads in each block fits into a np.uint32 and the total number 
    of basepairs fits into a np.uint64.
    
    Args:
        coverage_per_chrom: coverage per chromosome, as a dict of {chrom: np.ndarray of shape (L / blocksize,)}
        chrom_lens: chromosome lengths, as a dict of {chrom: int}
        blocksize: block size
    """
    def __init__(self, coverage_per_chrom, chrom_lens, blocksize):
        # note: given number of blocks, it is not possible to uniquely infer the chromosome length from the blocksize or vice versa (e.g. 4 blocks: chrom len = 35 -> blocksize = 9, 10, 11; blocksize = 10 -> chrom len = 31, 32,...)
        self.coverage_per_chrom = coverage_per_chrom
        assert all(_divide_round_up(chrom_len, blocksize) == len(self.coverage_per_chrom[chrom]) for (chrom, chrom_len) in chrom_lens.items())
        self.chrom_lens = chrom_lens
        if blocksize == 1:
            logger.warning("more efficient to use the class CoverageTracker")
        self.blocksize = blocksize
        
    @classmethod
    def empty_from_ref_genomes(cls, ref_genomes: Union[Dict[str, Union[os.PathLike, Number]], os.PathLike], blocksize=1, **kwargs):
        """
        Create empty coverage arrays (with shapes of the reference genomes) in the file
        
        For a dict {genome_name: genome_file}, it prefixes each chromosome name with the genome name.
        
        Args:
            ref_genomes: dict of {genome_name: genome_file} or path to a single genome file, or dict {chrom_name: chrom_len}
            blocksize: block size
            **kwargs: passed to constructor
        """
        coverage_per_chrom = {}
        chrom_lens = {}

        if isinstance(ref_genomes, dict):
            if isinstance(next(iter(ref_genomes.values())), Number):
                # provided lengths
                return cls.empty_from_lens(ref_genomes, blocksize=blocksize, **kwargs)
            
            # metagenome mode
            for genome_name, genome_path in ref_genomes.items():
                gen_coverage_per_chrom, gen_chrom_lens = cls._zero_arr_and_len_per_chrom(genome_path, seq_prefix=f"{genome_name}-", blocksize=blocksize)
                coverage_per_chrom.update(gen_coverage_per_chrom)
                chrom_lens.update(gen_chrom_lens)
        else:
            # genome mode: no seq prefix
            coverage_per_chrom, chrom_lens = cls._zero_arr_and_len_per_chrom(ref_genomes, blocksize=blocksize)
        
        return cls(coverage_per_chrom=coverage_per_chrom, chrom_lens=chrom_lens, blocksize=blocksize, **kwargs)
    
    ARRAY_DTYPE = np.uint32 # count > 255 is possible due to large blocksize
    @classmethod
    def empty_from_lens(cls, chrom_lens: Dict[str, int], blocksize=1, **kwargs):
        """
        Create empty coverage arrays according to lengths
        """
        coverage_per_chrom = {genome_name: np.zeros(_divide_round_up(genome_len, blocksize), dtype=cls.ARRAY_DTYPE) for (genome_name, genome_len) in chrom_lens.items()}
        return cls(coverage_per_chrom=coverage_per_chrom, chrom_lens=chrom_lens, blocksize=blocksize, **kwargs)
    
    @classmethod
    def _zero_arr_and_len_per_chrom(cls, genome_filename, seq_prefix="", blocksize=1) -> Dict[str, np.ndarray]:
        """
        Initialize an array with zeros for each chromosome in the genome
        
        Args:
            genome_filename: fasta file to identify chromosome names
            seq_prefix: prefix to add to each chromosome name
            blocksize: the block size to use for the coverage array, collapsing this many positions into one (except for the last)
            
        Returns:
            tuple of dict mapping chromosome to a zero array, dict mapping chromosome to its length
        """
        with pysam.FastaFile(genome_filename) as fasta:
            # use uint8 because counts should be less than 256 usually
            return (
                {seq_prefix + normalize_seq_name(id): np.zeros(_divide_round_up(fasta.get_reference_length(id), blocksize), dtype=cls.ARRAY_DTYPE) for id in fasta.references},
                {seq_prefix + normalize_seq_name(id): fasta.get_reference_length(id) for id in fasta.references}
            )

    def get_chrom_start_len(self, read_id) -> Optional[Tuple[Any, int, int]]:
        """
        Map a read to its chromosome, start and length on the reference
        
        Args:
            read_id: a NanoSim read id
        
        Returns:
            Tuple (chrom, ref_start, ref_len) with respect to forward strand; None
        """
        raise NotImplementedError()
    
    def add_reads(self, read_ids: Iterable, verbose=False):
        """
        Add new reads that were sequenced (as this influences coverage) 
        
        Args:
            read_ids: list or generator of read ids, should be parsable by NanoSimId.parse, 
                ref_len should be adapted if read was rejected in the middle
                
            verbose: whether to output regularly how many reads have been added
        
        Returns: 
            per chromosome, (number of reads added, number of basepairs added)
        """
        nb_reads = {chrom: 0 for chrom in self.coverage_per_chrom.keys()}
        nb_basepairs = {chrom: 0 for chrom in self.coverage_per_chrom.keys()}
        for read_id in tqdm.tqdm(read_ids) if verbose else read_ids:
            mapping_res = self.add_read(read_id)
            if mapping_res is not None:
                chrom, ref_start, ref_len = mapping_res
                nb_reads[chrom] += 1
                nb_basepairs[chrom] += ref_len
        
        if verbose:
            logger.info(f"Added {sum(nb_reads.values())} reads with {sum(nb_basepairs.values())} bps")
            
        return nb_reads, nb_basepairs
    
    def add_read(self, read_id):
        """
        Add a single read with this read_id.
        
        Circular genomes are not supported currently.
        """
        mapping_res = self.get_chrom_start_len(read_id)
        if mapping_res is None:
            logger.debug(f"Could not map read with id '{read_id}'")
            return None
        chrom, ref_start, ref_len = mapping_res
        # note: assumes that ref_len is read reference length, i.e. if read was rejected in the middle, the ref length was adapted accordingly
        # self.coverage_per_chrom[chrom][ref_start:ref_start+ref_len] += 1
        assert ref_start + ref_len <= self.chrom_lens[chrom], f"Read {read_id} mapped to {chrom}:{ref_start}-{ref_start+ref_len} but chromosome length is {self.chrom_lens[chrom]}, circular genome currently not supported"
        _add_count_to_blocked_array(self.coverage_per_chrom[chrom], ref_start, ref_start + ref_len, blocksize=self.blocksize)
        return chrom, ref_start, ref_len

    def _avg_cov_per_block(self, chrom) -> float:
        """
        Compute average coverage per basepair for each block in a chromosome, accounting for last block being shorter
        """
        return self.coverage_per_chrom[chrom] / self._block_sizes(chrom)
    
    def _block_sizes(self, chrom) -> np.ndarray:
        """Compute the sizes of each block in a chromosome"""
        res = np.ones_like(self.coverage_per_chrom[chrom]) * self.blocksize
        if self.chrom_lens[chrom] % self.blocksize != 0:
            res[-1] = self.chrom_lens[chrom] % self.blocksize
        return res
    
    def avg_chrom_coverage(self, chrom, metric="mean") -> float:
        """
        Compute average coverage of a chromosome
        
        For "mean", computes the average coverage of the basepairs. 
        For "median", it computes the median of the average coverage of basepairs per block (the last block is shorter)
        
        Args:
            chrom: chromosome
            metric: metric to use, either "mean" or "median"
            
        Returns:
            average coverage for the chromosome
        """
        if metric == "mean":
            return self.coverage_per_chrom[chrom].sum(dtype=np.uint64) / self.chrom_lens[chrom]
        else:
            return np.median(self._avg_cov_per_block(chrom))
        
    def get_fraction_cov_atleast(self, threshold, chroms: Optional[List]=None) -> float:
        """
        Compute the fraction of basepairs in blocks with average basepair coverage above a given threshold.
        
        Count the fraction of basepairs which are in blocks with average coverage above a given threshold, then normalize.
        
        Args:
            threshold: coverage threshold
            chroms: chromosomes to consider for fraction
        """
        if chroms is None:
            chroms = list(self.coverage_per_chrom.keys())
        if len(chroms) == 0:
            logger.warning("get_fraction_cov_atleast called with empty chroms, returning 1.0")
            return 1.0
        # last block can be shorter, so we also have to weight it differently
        return sum(((self._avg_cov_per_block(chrom) >= threshold) * self._block_sizes(chrom)).sum(dtype=np.uint64) for chrom in chroms) / sum(self.chrom_lens[chrom] for chrom in chroms)
        
    def get_chrom_lens(self) -> Dict[str, int]:
        """
        Get lengths of chromosomes
        """
        return self.chrom_lens
    
    def plot_state(self, plot_type, target_coverage=None, **kwargs):
        """
        Plot coverage per chromosome with one axis per chromosome, as a line plot (over chromosome position) or histogram
        
        Args:
            plot_type: "line" or "hist"
            target_coverage: if given, plot a dashed line at this coverage level
        """
        assert plot_type in ["line", "hist"]
        
        import matplotlib.pyplot as plt
        
        num_chroms = len(self.coverage_per_chrom)
        fig, axes = plt.subplots(ncols=num_chroms, figsize=(6, 1 + 3*num_chroms), squeeze=False)
        for (chrom, ax) in zip(self.coverage_per_chrom, axes.flat):
            covs_per_position = self._avg_cov_per_block(chrom)
            plot_fcn = functools.partial(plot_cov_per_position_lineplot, blocksize=self.blocksize) if plot_type == "line" else plot_covs_hist
            plot_fcn(covs_per_position, ax=ax, target_coverage=target_coverage, cov_mean=self.avg_chrom_coverage(chrom), cov_median=self.avg_chrom_coverage(chrom, metric="median"))
            if target_coverage is None:
                ax.set_title(chrom)
            else:
                fraction_at_least = self.get_fraction_cov_atleast(threshold=target_coverage, chroms=[chrom])
                ax.set_title(f"{chrom}: {fraction_at_least*100:.2f}%")
            
        [fig.delaxes(ax) for ax in axes.flat[num_chroms:]] # remove extra axes
        fig.suptitle("Coverage per basepair" if self.blocksize == 1 else "Average coverage per block")
        
        return fig

# todo2: mixin
# CovTrackerClass = CoverageTracker # todo2
CovTrackerClass = BlockwiseCoverageTracker # todo2

class NanoSimCoverageTracker(CovTrackerClass):
    """
    Track coverage by parsing the location from the NanoSim read ids
    
    NanoSim unaligned reads do not map.
    """
    def get_chrom_start_len(self, read_id):
        nanosim_id = NanoSimId.from_str(read_id)
        if nanosim_id.read_type == "unaligned":
            return None
        return (nanosim_id.chrom, nanosim_id.ref_pos, nanosim_id.ref_len)
    
class PafCoverageTracker(CovTrackerClass):
    """
    Track coverage by looking up the read ids in a PAF file
    """
    def __init__(self, paf_df, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.paf_df = paf_df
        
    @classmethod
    def empty_from_paf(cls, paf_file, **kwargs):
        """
        Create an empty PafCoverageTracker from a PAF file
        
        Reads in the PAF file and extracts reference genome lengths from the PAF file
        
        Args:
            paf_file: PAF file to read
        """
        logger.debug(f"Reading PAF file '{paf_file}'")
        paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=[0, 5, 6, 7, 8], names=["read_id", "target_name", "target_length", "target_start", "target_end"])
        logger.debug("Finished reading PAF file")
        
        logger.debug(f"""The PAF file mapped {sum(paf_df["target_name"] != "*") / len(paf_df):.3%} of the reads""")
        paf_df = paf_df[paf_df["target_name"] != "*"] # remove unmapped
        paf_df["target_name"] = paf_df["target_name"].astype("category")
        paf_df["target_length"] = paf_df["target_length"].astype(np.int64)
        paf_df["target_start"] = paf_df["target_start"].astype(np.int64)
        paf_df["target_end"] = paf_df["target_end"].astype(np.int64)
        
        chrom_lens = paf_df.groupby("target_name", observed=True)["target_length"].first().to_dict()
        logger.info(f"Identified chromosomes with following lengths from PAF: {chrom_lens}")
        
        assert paf_df["read_id"].is_unique, "read ids in PAF file must be unique, probably took PAF file from ReadFish which logs each time a chunk maps" # see https://github.com/lh3/minimap2/blob/master/FAQ.md to filter out chimeric alignments
        assert all(paf_df["target_end"] > paf_df["target_start"]), "circular genomes currently not supported"
        paf_df["nb_ref_bps"] = paf_df["target_end"] - paf_df["target_start"] # no +1 offset it seems
        paf_df.drop(["target_end", "target_length"], axis=1, inplace=True)
        paf_df.set_index("read_id", inplace=True)
        
        return cls.empty_from_lens(chrom_lens=chrom_lens, paf_df=paf_df, **kwargs)
    
    def get_chrom_start_len(self, read_id):
        if read_id not in self.paf_df.index:
            return None
        return self.paf_df.loc[read_id][["target_name", "target_start", "nb_ref_bps"]].values