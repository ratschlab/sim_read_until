"""
Manage simulation parameters such as bps_per_second, chunk_size, gap_samplers
"""

from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler
import numpy as np
from typing import Dict, Optional, Union

from simreaduntil.simulator.pore_model import PoreModel


class SimParams:
    """
    Encapsulates simulation parameters
        
    Args:
        gap_samplers: gap samplers for each channel
        bp_per_second: basepairs per second going through the pore (per channel)
        default_unblock_duration: extra delay to reject a read / unblock a pore, in seconds
        chunk_size: chunk size for selective sequencing ReadUntil (size of chunks when sending data; chunks are concatenated; last chunk has shorter size)
        seed: seed for random number generator, same state set on all channels
    """
    def __init__(self, gap_samplers: Dict[str, GapSampler], bp_per_second=450, default_unblock_duration=0.1, chunk_size=200, pore_model: Optional[PoreModel]=None, seed: Union[int, np.random.Generator]=0):
        self.set(gap_samplers=gap_samplers, bp_per_second=bp_per_second, default_unblock_duration=default_unblock_duration, chunk_size=chunk_size, seed=seed, pore_model=pore_model)

    def restrict_to_channels(self, channels, rand_state):
        """Subset SimParams to some channels"""
        return SimParams(gap_samplers={channel: self.gap_samplers[channel] for channel in channels}, bp_per_second=self.bp_per_second, default_unblock_duration=self.default_unblock_duration, chunk_size=self.chunk_size, seed=rand_state)
        
    def __repr__(self) -> str:
        # repr(random_state) is not very informative (does not show seed, so we store it separately and display it here)
        return f"""SimParams(bp_per_second={self.bp_per_second}, default_unblock_duration={self.default_unblock_duration}, chunk_size={self.chunk_size}, initial_seed={self._initial_seed}, n_channels={len(self.gap_samplers)})"""

    def set(self, *, gap_samplers: Dict[str, GapSampler]=None, bp_per_second=None, default_unblock_duration=None, chunk_size=None, pore_model=None, seed=None):
        """
        Set parameters, None values are ignored
        """
        if gap_samplers is not None:
            n_channels = len(gap_samplers)
            assert 1 <= n_channels, f"Need at least one channel, not {n_channels}"
            
            self.gap_samplers = gap_samplers
        if bp_per_second is not None:
            self.bp_per_second = bp_per_second
        if default_unblock_duration is not None:
            self.default_unblock_duration = default_unblock_duration
        if chunk_size is not None:
            self.chunk_size = chunk_size
        if pore_model is not None:
            self.pore_model = pore_model
        if seed is not None:
            if isinstance(seed, np.random.Generator):
                self.random_state = seed
                self._initial_seed = "unavailable"
            else:
                self._initial_seed = seed # don't use outside of this class, just for __repr__
                self.random_state = np.random.default_rng(seed)
        self._check_sim_params()
        
    @property
    def n_channels(self):
        """Number of channels"""
        return len(self.gap_samplers)

    def _check_sim_params(self):
        assert isinstance(self.bp_per_second, (int, float))
        assert self.bp_per_second > 0

        assert isinstance(self.chunk_size, int)
        assert self.chunk_size > 0

        assert isinstance(self.default_unblock_duration, (int, float))
        assert self.default_unblock_duration >= 0
        
        # assert isinstance(self.pore_model, PoreModel)

        assert isinstance(self.random_state, np.random.Generator)
        