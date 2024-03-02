"""
Simulate a collection of channels to synchronously start, stop, run mux scans, issue ReadUntil API calls on the channels
"""

#%%
import argparse
import datetime
import enum
import itertools
from pathlib import Path
import queue
from textwrap import dedent, indent
import threading
import time
import tqdm
from typing import Any, Dict, Generator, Optional, Tuple, Union, List
from threading import Thread
import logging
import warnings
import numpy as np
import pandas as pd    
import seaborn as sns
import matplotlib.pyplot as plt

from simreaduntil.shared_utils.plotting import FIGURE_EXT
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger
from simreaduntil.shared_utils.merge_axes import save_fig_and_pickle
from simreaduntil.shared_utils.nanosim_parsing import NanoSimId
from simreaduntil.shared_utils.plotting import make_tight_layout, rotate_xticks
from simreaduntil.shared_utils.utils import _cycle_list_deep
from simreaduntil.simulator.readpool import ReadPool, ReadPoolFromFile
from simreaduntil.simulator.readswriter import RotatingFileReadsWriter
from simreaduntil.simulator.channel import Channel, StoppedReceivingResponse, UnblockResponse, plot_channels
from simreaduntil.simulator.channel_stats import ChannelStats, combine_stats
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.simulator.utils import new_thread_name
from simreaduntil.shared_utils.timing import cur_ns_time
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler
from simreaduntil.shared_utils.utils import dill_dump, setup_logger_simple, tqdm_if_verbose, subset_dict
from simreaduntil.shared_utils.thread_helpers import ThreadWithResultsAndExceptions


logger = setup_logger_simple(__name__)
"""module logger"""

#%%
class InexistentChannelsException(Exception):
    """
    Trying to use channels that are not available
    
    Channels in the ReadUntil API are 1-based.
    """
    # define constructor to be able to pass message
    def __init__(self, message):
        super().__init__(message)
        
class ActionType(int, enum.Enum):
    """Type of action to perform"""
    Unblock = 0
    StopReceiving = 1
    
    def to_str(self):
        return {ActionType.Unblock: "unblock", ActionType.StopReceiving: "stop_receiving"}[self]

class ReadUntilDevice:
    """
    Base class to interact with a device that supports ReadUntil sequencing
    
    Functions:
    data: get_basecalled_read_chunks
    actions: unblock_read, stop_receiving_read
    simulation actions: start/stop acts on all channels, is_running
    
    These functions may be called in parallel (except for start/stop), so appropriate locking must be implemented.
    """
        
    @property
    def mk_run_dir(self) -> Union[Path, str]:
        """
        Directory where full reads are written to
        """
        raise NotImplementedError()
    
    @property
    def n_channels(self) -> int:
        """
        Number of channels
        """
        raise NotImplementedError()
    
    def _onebased_available_channels(self) -> List[int]:
        """
        Available channels to block/stop receiving, get chunks from, 1-based
        """
        return list(range(1, 1 + self.n_channels))
    
    def _check_channels_available(self, channel_subset: List[int]=None):
        """
        Check that the provided channels are available to block/stop receiving, get chunks from
        """
        if not set(channel_subset).issubset(self._onebased_available_channels()):
            raise InexistentChannelsException(message=f"Tried to use {channel_subset}, but only channels {self._onebased_available_channels()} are available (channels are 1-based!)")
    
    def start(self, **kwargs):
        """
        Start the sequencing.
        """
        raise NotImplementedError()
    
    def stop(self):
        """
        Stop the sequencing, e.g. write all partial reads to a file (flush).
        """
        raise NotImplementedError()
    
    def run_mux_scan(self, t_duration: float):
        """
        Starts a mux scan on all channels
        
        If a mux scan is already running, the duration can be updated with this method.
        
        Returns:
            number of reads rejected due to mux scan (0 when a mux scan is already running)
            
        Raises:
            SimulationNotRunningException: if simulation was not already running
        """
        raise NotImplementedError()
    
    @property
    def is_running(self) -> bool:
        """
        Whether the device is sequencing
        """
        raise NotImplementedError()
    
    def get_basecalled_read_chunks(self, batch_size=None, channel_subset=None) -> Generator[Any, None, None]:
        """
        Get available read chunks from the selected channels, from at most 'batch_size' channels
        
        Channels that are not running are ignored, so it also works when the simulation is over (and returns empty chunks).
        """
        raise NotImplementedError()
    
    def get_action_results(self, **kwargs) -> Generator[Tuple[Any, float, int, str, Any], Any, Any]:
        """
        Get new results of actions that were performed with unblock and stop_receiving (mux scans etc not included)
        
        Only the results that were not yet requested are returned.
        
        Returns:
            list of tuples (read_id, time, channel, action_type, result)
        """
        raise NotImplementedError()
    
    def unblock_read(self, read_channel, read_id, unblock_duration=None) -> Optional[bool]:
        """
        Unblock (reject) read on channel
        
        If the channel is not running, do not raise an exception, rather log it as missed.
        
        Args:
            read_channel: channel where read is
            read_id: read id to unblock
            unblock_duration: apply unblock voltage for a duration (default duration if not provided)
        
        Returns:
            whether the action was performed (not performed if the read was already over), or None (depending on the implementation)
            If None, use get_action_results() to get the result later.
        """
        raise NotImplementedError()
        
    def stop_receiving_read(self, read_channel, read_id) -> Optional[StoppedReceivingResponse]:
        """
        Stop receiving read on channel
        
        If the channel is not running, do not raise an exception, rather log it as missed.
        
        Args:
            read_channel: channel where read is
            read_id: read id to stop receiving
        
        Returns:
            result of the action, i.e., missed, already stopped or stopped, or None (depending on the implementation)
            If None, use get_action_results() to get the result later.
        """
        raise NotImplementedError()
    
    def unblock_read_batch(self, reads: List[Tuple], unblock_duration=None) -> List[Optional[bool]]:
        """
        Unblock a set of reads.

        Args:
            reads (tuple): List of (channel, read_id)
            unblock_duration (float): time in seconds to apply unblock voltage
        
        Returns:
            list with True/False whether action succeeded (if read was still active)
        """
        res = []
        for (channel, read_id) in reads:
            res.append(self.unblock_read(channel, read_id=read_id, unblock_duration=unblock_duration))
        return res
                
    def stop_receiving_batch(self, reads: List[Tuple]) -> List[Optional[StoppedReceivingResponse]]:
        """
        Stop receiving data from a set of reads

        Args:
            reads (tuple): List of (channel, read_id)

        Returns:
            list with action results
        """
        res = []
        for (channel, read_id) in reads:
            res.append(self.stop_receiving_read(channel, read_id=read_id))
        return res
    
    def send_message(self, msg: str, severity):
        """Send a message to the device"""
        logger.info(f"Device received message with severity {severity}: {msg}")
        
    def device_info(self, **kwargs) -> str:
        """
        Print information about the device, e.g. parameters, channel statuses etc."""
        raise NotImplementedError()
    
class SimulationRunningState(str, enum.Enum):
    """Current running state of the simulator"""
    Stopped = "SimStopped"
    Starting = "SimStarting" # between start and until thread to forward channels is started
    Running = "SimRunning" # when forwarding
    Stopping = "SimStopping" # to signal that forwarding thread should stop at the next occasion, called by .stop()
    
class ONTSimulator(ReadUntilDevice):
    """
    This simulates an ONT device with support for the ReadUntil API. 
    It sets up a number of channels, the channels fetch reads from the ReadPool and write finished reads to the ReadsWriter.
    
    The ReadUntil API is supported by the following functions:
    - get_basecalled_read_chunks: get new chunks from the channels for the currently active reads
    - unblock_read: reject a read
    - stop_receiving_read: stop receiving a read
    Note: The channels are 1-based in the ReadUntil API (though the class internally uses 0-based indexing).
    
    Forwarding the simulation means forwarding each channel from the last time it was forwarded to the specified time. This leads to a number of 
    reads being written.
    It supports two operating modes that should NOT be mixed:
    - async mode: run a thread to forward the simulation at regular time intervals, so reads get written once finished; 
    acceleration_factor can be used to speed up the simulation
    methods: .start(), .stop(), .is_running
    - sync mode: forward the simulation whenever it is called; blocks the main thread
    methods: .sync_forward(t), .sync_start(t), .sync_stop(t)
      
    The reads are not necessarily written in an ordered fashion because the channels are traversed in a fixed order, not by read finishing time (difficult because a 
    channel may write several reads at once when forwarded). If desired, implement this at the ReadsWriter level.
    
    The simulation parameters can be changed during the simulation. They will be automatically picked up.
    
    This class is thread-safe. start/stop can be called in parallel and acquire a lock. They return true/false whether they succeeded. The forward loop acquires
    the lock on each iteration when it is forwarding the channels, then releases it, so the simulation can be stopped.
    
    Note: The internal thread that forwards the simulation works because it only forwards channels which are thread-safe. When the simulation
    is stopped, the thread to forward the channels finishes and the channels are stopped then only.
    
    Arguments:
        read_pool: where reads come from
        reads_writer: reads writer, ideally of type RotatingFileReadsWriter with attribute .output_dir (used in .mk_run_dir attribute)
        sim_params: simulation parameters, can be modified during the simulation
        channel_status_filename: where to write combined channel status at regular intervals
        output_dir: output dir returned by self.mk_run_dir, this is where files will be put
    """
    def __init__(self, read_pool: ReadPool, reads_writer: RotatingFileReadsWriter, sim_params: SimParams, channel_status_filename: Optional[Union[str, Path]]=None, output_dir="<unavailable>"):
        logger.debug(f"Creating ONT device simulator")
        
        self._read_pool = read_pool
        self._reads_writer = reads_writer # shared across all channels
        self.sim_params = sim_params
        self._channels: List[Channel] = [Channel(channel_name, read_pool, reads_writer, sim_params=sim_params) for channel_name in sim_params.gap_samplers.keys()]
        
        self._channel_status_filename = channel_status_filename
        self._output_dir = output_dir
        
        # thread that forwards simulation at regular time intervals, may not be alive, so call .is_running() to check if simulation is currently running
        self._forward_sim_thread: Optional[ThreadWithResultsAndExceptions] = None
        
        # the sim state captures the current state of the simulation (started, stopped etc.). Before modifying the state, we get a lock, 
        # check we are in the right initial state and then transition to a new state
        self.lock_sim_state = threading.RLock() # lock to hold running state of simulation fixed
        self.sim_state = SimulationRunningState.Stopped
        
        self.action_queue = queue.Queue() # queue for actions to perform on the simulator to avoid lock contention, do it at the end of each forward
        self._action_results = []
        
    @property
    def mk_run_dir(self) -> Union[Path, str]:
        return self._output_dir
    
    def device_info(self, sim_params=True, channel_states=False) -> str:
        """
        Print information about the device
        """
        with self.lock_sim_state:
            res = f"Simulator over {self.n_channels} channels"
            if sim_params:
                res += "\n" + dedent(f"""\
                Sim parameters:
                    sim_params: {self.sim_params}
                    n_channels={self.n_channels}
                    read_pool={self._read_pool}
                    sim_state={self.sim_state}""")
            if channel_states:
                res += "\nChannel states:\n" + "\n".join(map(repr, self._channels))
            return res
    
    @property
    def n_channels(self) -> int:
        return len(self._channels)
    
    def __repr__(self) -> str:
        return self.device_info(sim_params=False, channel_states=True)
    
    def get_channel_stats(self, combined=False) -> Union[ChannelStats, List[ChannelStats], str]:
        """
        Get channel statuses, per channel
        
        This method is not thread-safe because it calls .stats on each channel, which is not thread-safe.
        If the method is called while the simulation is running, the channel statuses may be inconsistent.
        
        Args:
            combined: if True, combine stats over all channels
        Returns:
            list of channel statuses, or combined stats
        """
        with self.lock_sim_state:
            if combined:
                return combine_stats((channel.stats for channel in self._channels))
            else:
                return [channel.stats for channel in self._channels]
    
    def _forward_channels(self, t, delta=False, show_progress=False):
        """
        Forward all channels to time t, optionally by delta
        
        Args:
            t: time to forward to
            delta: if True, interpret t as time difference
            show_progress: show progress across channels with tqdm
        """
        for channel in tqdm_if_verbose(self._channels, verbose=show_progress, desc="Forwarding channels"): # doing this in parallel in threads is not useful due to Python GIL
            channel.forward(t, delta=delta)
            
    def _start_channels(self, t=None):
        """
        Start channels, i.e. place first elements in channel
        
        Args:
            t: start time, if None, start at current time
        """
        if t is None:
            t = cur_ns_time()
        for channel in self._channels:
            channel.start(t)
            
    def _stop_channels(self):
        """Stop all channels and flush any pending reads to disk etc."""
        for channel in self._channels:
            assert channel.is_running
            channel.stop()
        self._reads_writer.finish()
        
    def _all_channels_finished(self) -> bool:
        """Whether no reads are left and all channels have finished"""
        
        # premature check: if read pool still has reads, then not all channels have finished
        return self._read_pool.definitely_empty and all(channel.is_idle() for channel in self._channels)
            
    ########## asynchronous mode ##########
    
    def start(self, *args, **kwargs):
        """Start the simulation.
        
        This method starts a thread that forwards the simulation at regular intervals.
        The simulation stops if all channels are done (since no reads are left) or it is explicitly stopped
        The  execution is not deterministic because of the use of time.sleep().
        
        Args:
            args, kwargs: forwarded to _start
            
        Returns:
            Whether the simulation was started (i.e. it was not running)
        """
        
        with self.lock_sim_state:
            if self.sim_state != SimulationRunningState.Stopped:
                logger.warning(f"Simulation is already running, not starting it")
                return False
        
            self.sim_state = SimulationRunningState.Starting
            
            # don't set as daemon because reads in-progress need to be written to a file
            self._forward_sim_thread = ThreadWithResultsAndExceptions(
                target=self._forward_sim_loop, name=new_thread_name("simforw-{}"), args=args, kwargs=kwargs
            )
            logger.info("Starting the simulation")
            
            self._start_channels(t=0) # start channels here so that channels are definitely running once this function finishes (and not once the thread executes this)
            self._forward_sim_thread.start()
            
            return True
        
    def _forward_sim_loop(self, acceleration_factor=1.0, update_method="realtime", log_interval: int=10, stop_if_no_reads=True, **kwargs):
        """
        Helper method launched by .start() to forward the simulation.
        
        It forwards the simulation at regular time intervals.
        It always forwards the simulation by the same amount of time (provided the simulation parameters stay unchanged).
        It breaks once self._please_stop is set (at the next iteration) or all channels have finished (because no reads are left).
        The channel status (aggregated across all channels) is written to a file at regular intervals (if filename provided).
        
        For debugging, it can be invoked directly rather than through .start().
        
        Args:
            acceleration_factor: by how much to speed up the simulation, it is safe to do so as long as the simulation can keep up with it (warnings are printed if not)
            update_method: how to forward the simulation (taking into account acceleration)
                - "realtime": update simulation time so it matches real time (accounting for acceleration); relevant when a delay happens between iterations
                - "constant": always update simulation time by same amount
                
            log_interval: number of steps between logging elapsed simulator time the combined channel status, set to 0 for no logging
            stop_if_no_reads: whether to stop the simulation if no reads are left            
        """
        with self.lock_sim_state:
            assert self.sim_state == SimulationRunningState.Starting
            self.sim_state = SimulationRunningState.Running
            
        logger.debug("Simulator forward thread started...")
        logger.info(f"Device info: {self.device_info()}")
        assert acceleration_factor > 0, f"invalid acceleration_factor {acceleration_factor}"
        
        if self._channel_status_filename is not None:
            self._channel_status_fh = open(self._channel_status_filename, "w")
            print(ChannelStats.get_table_header(), file=self._channel_status_fh)
            
        def _log():
            logger.debug(f"Simulation has been running for real {cur_ns_time() - t_real_start:.2f} seconds with acceleration factor {acceleration_factor:.2f} (t_sim={t_sim:.2f}, i={i})")
                
            combined_channel_statuses = self.get_channel_stats(combined=True)
            logger.info("\nCombined channel status: " + str(combined_channel_statuses))
            
            if self._channel_status_filename is not None:
                print(combined_channel_statuses.get_table_line(), file=self._channel_status_fh)
        
        t_real_start = cur_ns_time()
        t_real_last_forward = t_real_start
        t_sim = 0
        i = 0
        _log()
        while (self.sim_state == SimulationRunningState.Running) and (not self._all_channels_finished()):
            i += 1
            
            # delay as necessary
            time_now = cur_ns_time()
            time_sleep = self._compute_delta_t_sim() / acceleration_factor + t_real_last_forward - time_now
            if time_sleep < 0:
                logger.warning(f"Simulation cannot keep up, delay: {-time_sleep} seconds (iteration {i}), desired delta_t between iterations: {self._compute_delta_t_sim() / acceleration_factor}s")
            else:
                logger.debug(f"Simulation sleeping for {time_sleep} seconds (iteration {i})")
                time.sleep(time_sleep)
                
            with self.lock_sim_state:
                assert self.sim_state != SimulationRunningState.Stopped
                    
                # get time only once the lock was acquired
                if update_method == "constant":
                    t_sim += self._compute_delta_t_sim()
                else:
                    assert update_method == "realtime"
                    t_sim = (cur_ns_time() - t_real_start) * acceleration_factor
                    
                logger.debug(f"Forwarding to time {t_sim}")
                t_real_last_forward = cur_ns_time()
                self._forward_channels(t_sim)
                self._process_actions()
            
            if (log_interval != -1) and (i % log_interval == 0):
                _log()
        
        logger.info(f"Simulation ended at time t_sim={t_sim:.2f}s, t_real={cur_ns_time() - t_real_start}, {i} iterations, delta_t_real={self._compute_delta_t_sim() / acceleration_factor}s, consumed {self._read_pool.nb_reads_returned} reads")
        _log()
        if self._channel_status_filename is not None:
            self._channel_status_fh.close()
        
        if self.sim_state == SimulationRunningState.Stopping:
            logger.info(f"Simforward thread: Ended because a stop request was received")
        else:
            logger.warning(f"Simforward thread: Ended at simulated time {t_sim:.2f} because no reads left")
            if stop_if_no_reads:
                self.stop(_join_thread=False)
            
    def _compute_delta_t_sim(self):
        """
        Compute time to advance simulation by, dynamically recomputed
        
        Returns:
            Length of one chunk in seconds (without acceleration)
        """
        return self.sim_params.min_chunk_size / self.sim_params.bp_per_second
    
    def stop(self, _join_thread=True):
        """
        Stop simulation by stopping thread and stopping all channels.
        
        It will stop at the next time step.
        
        Args:
            _join_thread: whether to join the thread, only for internal use (set to False if stopping from the thread itself)
            
        Returns:
            Whether the simulation was stopped (i.e. it was running and not in the process of being stopped)
            The simulation has not necessarily stopped when this method returns False, it only started the stopping process.
        """
        
        logger.info("Stop request received, stopping simulation...")
        
        with self.lock_sim_state:
            if self.sim_state != SimulationRunningState.Running:
                return False
            self.sim_state = SimulationRunningState.Stopping
            logger.info("Set simulation state to Stopping...")
        
        if _join_thread:
            self._forward_sim_thread.join()  # block, try hard for .cancel() on stream
            self._forward_sim_thread.raise_if_error()
            assert not self._forward_sim_thread.is_alive()
        self._forward_sim_thread = None
        
        with self.lock_sim_state:
            self._stop_channels()
            self.sim_state = SimulationRunningState.Stopped
            logger.info("Set simulation state to Stopped...")
            
        return True
        
    @property
    def is_running(self):
        """
        Whether the simulation is currently not stopped in asynchronous mode.
        
        This is a momentary state, it may change at any time.
        If the simulation was started in synchronous mode, this will always return False.
        
        Returns:
            Whether the simulation is currently not stopped
        """
        return self.sim_state != SimulationRunningState.Stopped
    
    ############## chunk related methods ##############
    
    def get_basecalled_read_chunks(self, batch_size=None, channel_subset=None) -> Generator[Tuple[Any], None, None]:
        """
        It permutes the channels and gets at most 'batch_size' from them.
        Channels with no new chunks are filtered out.
        
        If some channels are not running (because the simulation is stopped or is being stopped), the
        chunks from those channels are simply not returned, instead of raising an error
        
        Args:
            channel_subset: restrict to these channels (if provided), 1-based, influenced by batch_size
            batch_size: maximum number of channels to get reads from
        Yields:
            basecalled chunks from channels: (channel_nb+1, read_id, chunk, quality, estimated ref len of all chunks returned so far for this read)
            channel_nb+1 is 1-based (since real ONT devices are)
        """
        channel_subset = self._onebased_available_channels() if channel_subset is None else channel_subset
        self._check_channels_available(channel_subset)
        if batch_size is None:
            batch_size = len(channel_subset)
        
        nb_chunks = 0
        # permute channels so that all are equally likely to be chosen (due to batch size), not always the first ones
        for channel in self.sim_params.random_state.permutation(channel_subset):
            if nb_chunks >= batch_size:
                break
            chunks, read_id, estimated_ref_len_so_far = self._channels[channel-1].get_new_samples() # if simulation was already stopped (in between), just returns ""
            # ignore if no new chunks (e.g. if channel does not have a read currently)
            if len(chunks) > 0: # chunks is either str or array of raw signals
                nb_chunks += 1
                yield (channel, read_id, chunks, "noquality", estimated_ref_len_so_far)
                
    def get_raw_chunks(self, *args, **kwargs):
        """
        Get raw chunks from channels
        
        Pore model must be set in sim_params
        Currently not implemented via the gRPC API
        """
        for (channel, read_id, chunks, quality, estimated_ref_len_so_far) in self.get_basecalled_read_chunks(*args, **kwargs):
            yield (channel, read_id, self.sim_params.pore_model.to_raw(chunks), quality, estimated_ref_len_so_far)
                
    def get_action_results(self, clear=True) -> Generator[Tuple[Any, float, int, str, Any], Any, Any]:
        """
        Get action results of actions performed on simulator
        
        Not thread-safe if clear=True.
        
        Args:
            clear: whether to clear the action results list after getting it
            
        Returns:
            List of action results: (read_id, time, channel, action, result)
        """
        if clear:
            # clear action list
            res = self._action_results
            self._action_results = []
            return res
        else:
            return self._action_results.copy()
    
    def unblock_read(self, read_channel, read_id, unblock_duration=None):
        self.action_queue.put((ActionType.Unblock, (read_channel, read_id, unblock_duration)))
        
    def stop_receiving_read(self, read_channel, read_id):
        self.action_queue.put((ActionType.StopReceiving, (read_channel, read_id)))
        
    # process actions asynchronously after calling forward since otherwise, there is a lot of lock contention which
    # means we cannot run at acceleration factor 10
    def _process_actions(self):
        while True:
            try:
                action_type, args = self.action_queue.get_nowait()
            except queue.Empty:
                return
            if action_type == ActionType.Unblock:
                self._unblock_read(*args)
            else:
                assert action_type == ActionType.StopReceiving
                self._stop_receiving_read(*args)
        
    def _unblock_read(self, read_channel, read_id, unblock_duration=None) -> Optional[bool]:
        """Unblock read"""
        self._check_channels_available([read_channel])
        action_res = self._channels[read_channel-1].unblock(unblock_duration=unblock_duration, read_id=read_id)
        self._action_results.append((read_id, self._channels[read_channel-1].t, read_channel, ActionType.Unblock, action_res))
        # logger.info(f"Unblocking read {read_id} on channel {read_channel}, result: {action_res.to_str()}")
        return action_res
        
    def _stop_receiving_read(self, read_channel, read_id) -> Optional[StoppedReceivingResponse]:
        """Stop receiving from read"""
        self._check_channels_available([read_channel])
        action_res = self._channels[read_channel-1].stop_receiving(read_id=read_id)
        self._action_results.append((read_id, self._channels[read_channel-1].t, read_channel, ActionType.StopReceiving, action_res))
        # logger.info(f"Stopping receiving from read {read_id} on channel {read_channel}, result: {action_res.to_str()}")
        return action_res
    
    def run_mux_scan(self, t_duration: float, is_sync=False) -> int:
        """Pass in duration on each channel rather than end time because the channel may already have been forwarded in-between"""
        with self.lock_sim_state:
            # the lock ensures that channels are not forwarded
            if (not is_sync) and self.sim_state != SimulationRunningState.Running:
                logger.warning("Simulation not (or no longer) running, mux scan ignored")
                return 0
            if self._channels[0].has_active_mux_scan():
                # only check first channel because mux scans are almost synced
                logger.warning("mux scan already running, adding extra time")
            logger.info(f"Starting mux scan at sim time {self._channels[0].t:.2f}s")
            nb_rejected_reads = 0
            for channel in self._channels:
                nb_rejected_reads += int(channel.run_mux_scan(t_duration=t_duration))
            return nb_rejected_reads
            
    ############## synchronous mode ##############
    
    # helps to ensure that synchronous and asynchronous mode are not mixed
    def _check_not_async_mode(self):
        assert not self.lock_sim_state == SimulationRunningState.Stopped, "should not combine both methods (synchronous + asynchronous forward with thread)"
        
    def sync_forward(self, t, delta=False, show_progress=False):
        """
        Process actions and forward all channels to time t
        
        Using (t=0, delta=True) means actions are processed
        """
        self._check_not_async_mode()
        self._process_actions()
        return self._forward_channels(t, delta=delta, show_progress=show_progress)
    
    def sync_start(self, t=None):
        """
        Start all channels at time t
        
        Args:
            t: if None, starts at current time (slightly differs between channels as they are sequentially started)
        """
        self._check_not_async_mode()
        self._start_channels(t=t)
        
    def sync_stop(self):
        """
        Stop all channels
        """
        self._check_not_async_mode()
        self._stop_channels()
    
    ############## plotting channels ##############
    
    @property
    def save_elems(self):
        """Whether channels are currently saving their history (needed for plotting)"""
        return self._channels[0].save_elems
    @save_elems.setter
    def save_elems(self, val):
        """Set all channels to save their history (needed for plotting)"""
        for channel in self._channels:
            channel.save_elems = val
    
    def plot_channels(self, channel_subset=None, *args, **kwargs):
        """
        Plot (a subset of) channels
        
        Only those elements that were happening while .save_elems was True are plotted.
        
        Args:
            channel_subset: 1-based, if None, all channels are plotted
        """
        if channel_subset is None:
            channel_subset = self._onebased_available_channels()
        else:
            self._check_channels_available(channel_subset)
        return plot_channels([self._channels[idx-1] for idx in channel_subset], *args, **kwargs)
    
    ############## support pickling ##############
    def __getstate__(self):
        state = self.__dict__.copy()
        del state["lock_sim_state"]
        return state
    
    def __setstate__(self, state):
        # must properly set the _read_pool and _reads_writer after pickling
        self.__dict__.update(state)
        self.lock_sim_state = threading.RLock()
        
def convert_action_results_to_df(action_results):
        """
        Convert action results returned by simulator to a dataframe
        
        Args:
            action_results: returned by get_action_results
        
        Returns:
            data frame    
        """
        action_results = [(read_id, time, channel, action_type.to_str(), bool(action_res)) for (read_id, time, channel, action_type, action_res) in action_results]
        # (read_channel, read_id, ActionType.Unblock, action_res)
        action_results_df = pd.DataFrame.from_records(action_results, columns=["read_id", "time", "channel", "action_type", "success"])
        return action_results_df
    
def write_simulator_stats(simulators, output_dir=None):
    """
    Dump action results (success or missed) and channel statistics to the run dir
    
    If the action results were previously queried with clear=True, this list is incomplete.
    
    Args:
        simulators: simulators to combine (action results and stats are combined)
        output_dir: output directory; if None, same directory as reads
    """
    if len(simulators) == 0:
        return
    
    if output_dir is None:
        output_dir = simulators[0].mk_run_dir
    
    action_results_df = pd.concat([convert_action_results_to_df(simulator.get_action_results(clear=False)) for simulator in simulators], ignore_index=True)
    action_results_df.to_csv(output_dir / "action_results.csv", sep="\t", index=False)
    # plot_sim_actions(action_results_df)
    
    dill_dump(list(itertools.chain.from_iterable(simulator.get_channel_stats() for simulator in simulators)), output_dir / "channel_stats.dill")
    
def plot_sim_actions(action_results_df, save_dir=None, close_figures=None):
    """
    Plot simulator actions statistics
    """
    if close_figures is None:
        close_figures = save_dir is not None
        
    def close_fig(fig):
        if close_figures:
            plt.close(fig)
         
    fig = plot_nb_actions_per_read(action_results_df); logger.debug("Created 1 plot"); close_fig(fig)
    fig = plot_action_success_rate(action_results_df); logger.debug("Created 1 plot"); close_fig(fig)
            
def plot_nb_actions_per_read(action_results_df, save_dir=None):
    """
    Plot number of actions, number of unique per read
    
    Also issues a logger warning for those reads with multiple actions or contradicting actions.
    
    Args:
        action_results_df: dataframe with columns read_id, action_type
        save_dir: if not None, save figures to this directory
    """
    nb_actions_per_read = action_results_df.groupby("read_id", observed=True).size()
    # print(nb_actions_per_read.value_counts())

    nb_unique_actions_per_read = action_results_df.groupby("read_id", observed=True)["action_type"].nunique()
    # print(nb_unique_actions_per_read.value_counts())

    # by action type, count number of actions per read
    nb_actions_per_read_and_type = action_results_df.groupby(["read_id", "action_type"], observed=True).size()
    # print(nb_actions_per_read_and_type.groupby("action_type").value_counts(), observed=True)


    reads_with_multiple_actions = nb_actions_per_read[nb_actions_per_read > 1].index.values
    reads_with_contradicting_actions = nb_unique_actions_per_read[nb_unique_actions_per_read > 1].index.values
    if len(reads_with_multiple_actions) > 0:
        logger.warning(f"There are {len(reads_with_multiple_actions)} reads with multiple actions (possible same actions): {reads_with_multiple_actions}")
    if len(reads_with_contradicting_actions) > 0:
        logger.warning(f"There are {len(reads_with_contradicting_actions)} reads with contradicting actions (e.g. stop_receiving and unblock): {reads_with_contradicting_actions}")


    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 4))
    nb_actions_per_read.value_counts().plot.bar(ax=ax1)
    # nb_actions_per_read.plot.hist(ax=ax1)
    # sns.barplot(nb_actions_per_read.value_counts().to_frame().reset_index(names="nb_actions"), x="nb_actions", y="count", ax=ax1)
    ax1.set_yscale("log")
    ax1.set_xlabel("Number of actions")
    ax1.set_ylabel("Number of reads")
    ax1.set_title("Number of actions per read")

    nb_unique_actions_per_read.value_counts().plot.bar(ax=ax2)
    ax2.set_yscale("log")
    ax2.set_xlabel("Number of unique actions")
    ax2.set_ylabel("Number of reads")
    ax2.set_title("Number of unique actions per read")

    sns.barplot(nb_actions_per_read_and_type.groupby("action_type", observed=True).value_counts().to_frame().reset_index(names=["action_type", "nb_actions"]), x="nb_actions", y="count", hue="action_type", ax=ax3)
    ax3.set_yscale("log")
    ax3.set_xlabel("Number of actions")
    ax3.set_ylabel("Number of reads")
    ax3.set_title("Number of actions per read and type")
    
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"nb_actions_per_type.{FIGURE_EXT}")
    
    return fig

def plot_action_success_rate(action_results_df, save_dir=None):
    """Plots action success rate over time"""
    assert len(action_results_df) > 0, "no actions performed"
    
    action_results_df = action_results_df.sort_values("time")
    
    action_results_df["cum_nb_success_per_type"] = action_results_df.groupby("action_type", observed=True)["success"].cumsum()
    action_results_df["cum_nb_actions_per_type"] = action_results_df.groupby("action_type", observed=True).cumcount() + 1
    action_results_df["cum_success_rate_per_type"] = action_results_df["cum_nb_success_per_type"] / action_results_df["cum_nb_actions_per_type"]
    
    action_results_df = action_results_df.sample(n=min(len(action_results_df), 500))
    
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 5))
    sns.lineplot(data=action_results_df, x="cum_nb_actions_per_type", y="cum_nb_success_per_type", hue="action_type", ax=ax1)
    ax1.axline((0, 0), slope=1, color="black", linestyle="--", label="x=x")
    ax1.set_xlabel("Number of actions")
    ax1.set_ylabel("Number of successful actions")
    
    sns.lineplot(data=action_results_df, x="time", y="cum_nb_success_per_type", hue="action_type", ax=ax2)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Number of successful actions")
    rotate_xticks(ax2, rotation=45)
    
    sns.lineplot(data=action_results_df, x="time", y="cum_success_rate_per_type", hue="action_type", ax=ax3)
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Fraction of successful actions")
    rotate_xticks(ax3, rotation=45)
    
    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"action_success_rate.{FIGURE_EXT}")
    
    return fig

# class ReadUntilClientFromDevice(ReadUntilDevice):
#     """
#     ReadUntilClient with ReadUntil actions operating on a batch of reads
    
#     Named ReadUntilClientFromDevice to avoid nameclash with ReadUntilClient
    
#     start, stop, device_info not implemented. They should be directly called on the device.
#     """
#     def __init__(self, device : ReadUntilDevice):
#         self._device = device
        
#     def __repr__(self):
#         res = "ReadUntilClientFromDevice of the following device:\n"
#         res += indent(repr(self._device), "    ")
#         return res
    
#     @property
#     def n_channels(self) -> int:
#         """
#         Number of channels
#         """
#         return len(self._device.n_channels)
    
#     @property
#     def is_running(self) -> bool:
#         """
#         Whether the device is sequencing
#         """
#         return self._device.is_running
    
#     @property
#     def mk_run_dir(self):
#         return self._device.mk_run_dir
    
#     def get_basecalled_read_chunks(self, batch_size=None, channel_subset=None):
#         """
#         Yield basecalled chunks from channels
        
#         Args:
#             batch_size: maximum number of channels to get reads from
#             channel_subset: restrict to these channels (if provided)
        
#         Yields:
#             basecalled chunks from channels in the form (chan_key, read_id, chunk, quality, estimated ref len of all chunks returned so far for this read)
#         """
#         yield from self._device.get_basecalled_read_chunks(batch_size, channel_subset)
        
#     def unblock_read(self, read_channel, read_id, unblock_duration=None) -> bool:
#        return self._device.unblock_read(read_channel, read_id=read_id, unblock_duration=unblock_duration)
        
#     def stop_receiving_read(self, read_channel, read_id) -> StoppedReceivingResponse:
#         return self._device.stop_receiving_read(read_channel, read_id=read_id)
    
def stop_simulation_after_time_thread(simulator: ONTSimulator, t: float):
    """
    Stop simulation after some time
    
    Note: This may result in a race condition (.stop()) if the time limit is reached and the user simultaneously interrupts the simulation with Ctrl+C
    
    Args:
        simulator: simulator to stop
        t: time in seconds, real time (not sim time)
        
    Returns:
        thread to start, must be started with .start()!
    """
    def _stop_after_time():
        logger.info(f"Simulation will be stopped after {t:.2f}s (real time)")
        time.sleep(t)
        logger.info(f"Stopping simulation due to time limit of {t:.2f}s (real time)")
        return simulator.stop()
        
    x = ThreadWithResultsAndExceptions(target=_stop_after_time, daemon=True, name="stop_sim_after_time")
    return x

def run_periodic_mux_scan_thread(simulator: ONTSimulator, period: float, scan_duration: float, acceleration_factor: float=1):
    """
    Run periodic mux scan at times [period, 2*period, ...] each time for a given duration
    
    It stops once the simulation is no longer running.
    
    Args:
        simulator: simulator to run mux scan on
        period: how often to run mux scan, in seconds, sim time (not real time)
        scan_duration: mux scan duration, in seconds, sim time (not real time)
        acceleration_factor: acceleration factor
        
    Returns:
        thread to start, must be started with .start()!
    """
    assert period > 0
    assert scan_duration >= 0
    
    if scan_duration > period - 0.05:
        warnings.warn(f"Period between mux scans may be so short that mux scans happen the whole time: scan_duration={scan_duration:.2f}s, period={period:.2f}s")
    
    def _run_periodic_mux_scan():
        logger.info(f"Running periodic mux scan every {period:.2f}s (sim time) with duration {scan_duration:.2f}s, acceleration factor {acceleration_factor:.2f}")
        i = 1
        time_start = cur_ns_time()
        while True:
            time.sleep(period / acceleration_factor)
            if simulator.is_running:
                logger.info(f"Starting {i}th mux scan for {scan_duration:.2f}s (sim time), real time since start: {cur_ns_time() - time_start:.2f}s")
                simulator.run_mux_scan(t_duration=scan_duration)
            else:
                break
            i += 1
        logger.info("Stopped periodic mux scan because the device is no longer running")
        
    x = ThreadWithResultsAndExceptions(target=_run_periodic_mux_scan, daemon=True, name="periodic_mux_scan")
    return x

# todo2: move this function to the other read generation functions
from simreaduntil.shared_utils.dna import get_random_DNA_seq
def _random_seqs_gen_from_lens_with_cycle(lengths, channel, cycle_read_durations, use_nanosim_id, random_state=np.random.default_rng(2)):
    # cycling if read durations are exhausted
    cycle_fcn = itertools.cycle if cycle_read_durations else lambda x: x
    for (i, read_len) in enumerate(cycle_fcn(itertools.chain(lengths, ["end"]))):
        if read_len == "end":
            logger.info(f"Cycle of read durations for channel {channel} exhausted" + (", starting new cycle" if cycle_read_durations else ""))
            continue
        yield (
            str(NanoSimId(chrom="chr1", ref_pos=1, read_nb=f"read-ch{channel}-{i+1}", direction="F", ref_len=read_len)) 
            if use_nanosim_id else f"read_ch{channel}_{i+1}", 
            get_random_DNA_seq(n=int(read_len), random_state=random_state)
        )
        
def run_simulator_from_sampler_per_channel(
    mk_run_dir, sim_params: SimParams, 
    read_durations_per_channel: Dict[Any, Generator[Any, None, None]], cycle_read_durations: bool,
    seq_end_time=None, random_state=np.random.default_rng(2), use_nanosim_id=False
):
    """
    Run the ONTSimulator using gaps from a gap sampler, until maximum end time of the gap samplers.
    
    The simulation is replicated channelwise.
    Note: This means that the reads are not ordered by time.
    
    The read durations are reccycled if exhausted.

    Args:
        mk_run_dir: where to put the run outputs (reads)
        sim_params: simulator parameters
        read_durations_per_channel: durations of reads per channel
        cycle_read_durations: whether to cycle read durations if exhausted
        seq_end_time: how long to run the simulation; if None, use the maximum end time of the gap samplers
        random_state: random state
        use_nanosim_id: whether to use a NanoSim id (with ref_len = seq_len) rather than simple incrementing integers
        
    Returns:
        stopped simulator, filename to reads FASTA
    """
    from simreaduntil.simulator.readpool import ReadPoolFromIterablePerChannel
    from simreaduntil.simulator.readswriter import SingleFileReadsWriter
    from simreaduntil.shared_utils.dna import get_random_DNA_seq

    assert set(sim_params.gap_samplers.keys()) == set(read_durations_per_channel.keys()), "keys must agree"

    if seq_end_time is None:
        seq_end_time = max(gap_sampler.full_duration for gap_sampler in sim_params.gap_samplers.values())
    logger.info(f"Simulating for time {seq_end_time}s")

    # not fully exact because reads have different speeds in reality, pick the number of basepairs to best match the read duration
    read_pool = ReadPoolFromIterablePerChannel({
        channel: _random_seqs_gen_from_lens_with_cycle(
            np.rint(read_durations * sim_params.bp_per_second), channel=channel, cycle_read_durations=cycle_read_durations, 
            use_nanosim_id=use_nanosim_id, random_state=random_state
        ) for (channel, read_durations) in read_durations_per_channel.items()
    })

    # note: can be parallelized by splitting channels into subgroups and then use an individual simulator, reads_writer and 
    # read_pool per subgroup, also set the random state so that it is not always the same; however, this speed improvement
    mk_run_dir.mkdir(exist_ok=True)
    reads_filename = mk_run_dir / "reads.fasta"
    with open(reads_filename, "w") as reads_file:
        reads_writer = SingleFileReadsWriter(reads_file)

        simulator = ONTSimulator(
            read_pool=read_pool,
            reads_writer=reads_writer,
            sim_params=sim_params,
            output_dir=reads_writer.output_dir,
        )

        simulator.sync_start(0)
        # making one huge forward is fine because read durations are per channel, so there is no difference in how reads are 
        # chosen depending on whether all channels are forwarded by a small amount many times vs each channel once by a large time
        # going until original end time should always work because the replicated run is not shorter than the original run
        simulator.sync_forward(seq_end_time, show_progress=True)
        simulator.sync_stop()
        
    return [(simulator, reads_filename)]

def _get_spawnable_random_state(seed):
    # also support SeedSequence (not just Generator) to support py38 and np>=1.22 (np1.25 requires py39) where Generator.spawn does not exist
    res = np.random.default_rng(seed)
    if hasattr(res, "spawn"):
        return res
    return np.random.SeedSequence(seed)
    
def run_simulator_from_sampler_per_channel_parallel(
    mk_run_dir, sim_params: SimParams, 
    read_durations_per_channel: Dict[Any, Generator[Any, None, None]], cycle_read_durations: bool,
    seq_end_time=None, random_state: Union[np.random.Generator, np.random.SeedSequence]=_get_spawnable_random_state(2), use_nanosim_id=False, parallel_args=None
):
    """
    Run the ONTSimulator using gaps from a gap sampler in parallel, until maximum end time of the gap samplers.
    
    The simulation is replicated channelwise.
    Note: This means that the reads are not ordered by time.
    
    The read durations are reccycled if exhausted.
    
    Performance on my Mac: For 50 channels for 230_000s, time to forward all channels:
    - running joblib with 1 process: 03m06s
    - running joblib with -1 (12) processes: 33s
    - running without joblib as 1 process (run_simulator_from_sampler_per_channel_nonparallel): 02m57s

    Args:
        mk_run_dir: where to put the run outputs (reads)
        sim_params: simulator parameters
        read_durations_per_channel: durations of reads per channel
        cycle_read_durations: whether to cycle read durations if exhausted; prints a message whenever a new cycle starts
        seq_end_time: how long to run the simulation; if None, use the maximum end time of the gap samplers
        random_state: random state
        use_nanosim_id: whether to use a NanoSim id (with ref_len = seq_len) rather than simple incrementing integers
        parallel_args: arguments to pass to joblib.Parallel; if None or '"n_jobs": 1', joblib does not run in parallel (default if not set)
        
    Returns:
        stopped simulator, filename to reads FASTA
    """
    from simreaduntil.simulator.readpool import ReadPoolFromIterablePerChannel
    from simreaduntil.simulator.readswriter import SingleFileReadsWriter
        
    root_logger_level = logging.getLogger().level
    
    def forward_channels(idx, sim_params, read_durations_per_channel, random_state):
        # temporary hack because logger not inherited into new process (e.g. with loky backend) due to global state of loggers, see https://github.com/joblib/joblib/issues/1017
        if not logging.getLogger().hasHandlers():
            add_comprehensive_stream_handler_to_logger(None)
            logging.getLogger().setLevel(root_logger_level)
            # copying loggers will not register it with the logging framework, so we need to recreate them to set their level
        logger = setup_logger_simple(__name__)
        
        # forward a set of channels
        assert set(sim_params.gap_samplers.keys()) == set(read_durations_per_channel.keys()), "keys must agree"
        
        # logger.info(f"Starting simulation for idx {idx}")
        # not fully exact because reads have different speeds in reality, pick the number of basepairs to best match the read duration
        read_pool = ReadPoolFromIterablePerChannel({
            channel: _random_seqs_gen_from_lens_with_cycle(
                np.rint(read_durations * sim_params.bp_per_second), channel=channel, cycle_read_durations=cycle_read_durations, 
                use_nanosim_id=use_nanosim_id, random_state=random_state
            ) for (channel, read_durations) in read_durations_per_channel.items()
        })

        reads_filename = mk_run_dir / f"reads{idx}.fasta"
        with open(reads_filename, "w") as reads_file:
            reads_writer = SingleFileReadsWriter(reads_file)

            simulator = ONTSimulator(
                read_pool=read_pool,
                reads_writer=reads_writer,
                sim_params=sim_params,
                output_dir="<unavailable>",
            )

            simulator.sync_start(0)
            # forwarding channels to the end one at a time works because read durations are per channel, so there is no 
            # difference in how reads are chosen (compared to when channels are moved by a small amount many times)
            simulator.sync_forward(seq_end_time)
            simulator.sync_stop()
            return simulator, reads_filename
    
    if seq_end_time is None:
        seq_end_time = max(gap_sampler.full_duration for gap_sampler in sim_params.gap_samplers.values())
    logger.info(f"Simulating for time {seq_end_time}s")
    
    mk_run_dir.mkdir(exist_ok=True)
    assert set(sim_params.gap_samplers.keys()) == set(read_durations_per_channel.keys()), "keys must agree"
    channels = list(read_durations_per_channel.keys())
    
    spawned_random_states = random_state.spawn(len(channels))
    # also support SeedSequence (not just Generator) to support py38 and np=1.22 (np1.25 requires py39) where generator.spawn does not exist
    if isinstance(random_state, np.random.SeedSequence):
        spawned_random_states = [np.random.default_rng(rs) for rs in spawned_random_states]
    del random_state
    
    # parallelize by using an individual simulator, reads_writer and read_pool for each channel, 
    # also set the random state, otherwise may get copied
    parallel_args = {} if parallel_args is None else parallel_args
    from joblib import Parallel, delayed, effective_n_jobs
    logger.info(f"Running joblib with {effective_n_jobs(parallel_args.get('n_jobs', None))} workers")
    from joblib.externals.loky import set_loky_pickler
    set_loky_pickler("dill")
    # set_loky_pickler("cloudpickle")
    with Parallel(**parallel_args) as parallel:
        res = parallel(delayed(forward_channels)(
            idx=channel, sim_params=sim_params.restrict_to_channels([channel], rand_state=child_state), 
            read_durations_per_channel=subset_dict(read_durations_per_channel, [channel]),
            random_state=child_state,
        ) for (channel, child_state) in zip(tqdm.tqdm(channels, desc="Simulating channels"), spawned_random_states))
        
    return res

def get_simulator_delay_over_time_df(log_filename):
    """Extract the delay of the simulator over time from the log file"""
    MARKER = "Simulation cannot keep up, delay: "
    def parse_line(line):
        # return time of log entry, number of reads, time elapsed
        log_time, remaining = line.split(" - ", maxsplit=1)
        # convert log_time to time
        log_time = datetime.datetime.strptime(log_time, "%Y-%m-%d %H:%M:%S,%f")
        remaining = remaining.split(" --- ", maxsplit=1)[0]
        remaining = remaining.split(MARKER, maxsplit=1)[1]
        delay, remaining = remaining.split(" seconds (iteration ")
        iter_nb = remaining.split("), desired")[0]
        return log_time, float(delay), int(iter_nb)

    # line = "2023-08-21 15:10:55,471 - Simulation cannot keep up, delay: 0.06818270683288574 seconds (iteration 12), desired delta_t between iterations: 0.05400137741046832s --- simulator.py:464 (_forward_sim_loop) WARNING ##"
    # parse_line(line)
    
    with open(log_filename) as f:
        info = [parse_line(line) for line in f if MARKER in line]
        # info = list(itertools.islice((parse_line(line) for line in f if MARKER in line), 100))
    df = pd.DataFrame.from_records(info, columns=["time", "delay", "iteration"])
    if len(df) > 0:
        df["time"] = (df["time"] - df["time"].iloc[0]).dt.total_seconds()
    
    return df

def plot_simulator_delay_over_time(df, n_delays=200, save_dir=None):
    """Plot simulator delay over time for the loop updating the channels"""

    # df = df.sample(min(len(df), 200))
    # restrict to the largest rather than sample
    df = df.sort_values("delay", ascending=False).iloc[:n_delays]

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(11, 4))
    
    sns.lineplot(df, x="time", y="delay", ax=ax1)
    ax1.set_xlabel("Real time (of iteration) (s)")
    ax1.set_ylabel("Delay (s)")
    sns.lineplot(df, x="iteration", y="delay", ax=ax2)
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Delay (s)")
    
    fig.suptitle(f"Simulator loop delays (largest {n_delays})")

    make_tight_layout(fig)
    if save_dir is not None:
        save_fig_and_pickle(fig, save_dir / f"simulator_delay.{FIGURE_EXT}")

    return fig

def assign_read_durations_to_channels(read_durations_per_channel: List[np.ndarray], channel_names) -> Dict[Any, np.ndarray]:
    """
    Attribute read durations to channels
    
    If there are more channels than read_durations entries, the list is recycled (copy).
    Note: Cannot just return empty durations (since a channel may have reads rather than being blocked). Read duration
    of 0 is also not allowed since the sequence must have length > 0
    
    Args:
        read_durations_per_channel: list with each entry being a list of read durations of a channel
        channel_names: channel names, read durations will be assigned to them in this order
        
    Returns:
        dictionary {channel -> read durations}
    """
    # res = {chan: np.array([fallback_read_duration]) for chan in channel_names}
    # res.update(
    #     {chan: read_durations for (chan, read_durations) in zip(channel_names, read_durations_per_channel)}
    # )
    # return res
    return {chan: read_durations for (chan, read_durations) in zip(channel_names, _cycle_list_deep(read_durations_per_channel))}
    