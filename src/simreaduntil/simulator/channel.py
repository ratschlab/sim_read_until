"""
Model a single channel as a finite-state machine of channel elements
"""

#%%
from __future__ import annotations
import contextlib
import enum
from threading import Lock # for referring in type hints of a class's method to class itself, Python 3.7, otherwise use strings; see https://stackoverflow.com/questions/55320236/does-python-evaluate-type-hinting-of-a-forward-reference

from typing import Iterable, List, Union, Tuple, Any
from matplotlib.layout_engine import TightLayoutEngine
import numpy as np
from simreaduntil.shared_utils.logging_utils import setup_logger_simple
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler
from simreaduntil.simulator.simulator_params import SimParams

from simreaduntil.simulator.utils import in_interval
from simreaduntil.simulator.readpool import NoReadLeftException, ReadPool
from simreaduntil.simulator.readswriter import ReadsWriter
from simreaduntil.simulator.channel_element import ChannelBroken, ChannelElement, ShortGap, MuxScan, NoReadLeftGap, UnblockDelay, ChunkedRead, LongGap, ReadEndReason

from simreaduntil.shared_utils.plotting import get_fig_ax
#%%

__all__ = ["Channel", "ChannelNotRunningException", "ChannelAlreadyRunningException"]

logger = setup_logger_simple(__name__)
"""module logger"""

DONT_WRITE_ZERO_LENGTH_READS = True
"""Whether to pass finished reads with length 0 to the reads writer or not"""

class ChannelNotRunningException(Exception):
    """
    When the channel is not running and something requires the channel to be running
    """
    pass

class ChannelAlreadyRunningException(Exception):
    """
    When the channel is already running and something requires the channel to not be running
    """
    pass
    
class StoppedReceivingResponse(int, enum.Enum):
    """
    Action response when stopping to receive a read
    
    Can be used as a boolean to check whether the read is set to stopped receiving or not
    Convertible to boolean
    """
    MISSED = 0 # no longer current read
    STOPPED_RECEIVING = 1 # stopped receiving, not already in this state before
    ALREADY_STOPPED_RECEIVING = 2 # already stopped receiving
    
    def to_str(self):
        return {StoppedReceivingResponse.MISSED: "missed", StoppedReceivingResponse.STOPPED_RECEIVING: "stopped_receiving", StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING: "already_stopped_receiving"}[self]
    
class UnblockResponse(int, enum.Enum):
    """
    Action response of an unblock, convertible to boolean
    """
    MISSED = 0 # no longer current read
    UNBLOCKED = 1 # successfully unblocked
    
    def to_str(self):
        return {UnblockResponse.MISSED: "missed", UnblockResponse.UNBLOCKED: "unblocked"}[self]
    
class Channel:
    """
    Simulate the reads from a flow cell pore (channel)
    
    It keeps track of the current element in the channel, either a read, a 
    gap between reads, an unblocking delay, a channel blockage,
    no read left gap (infinite length).
    A read is passed to the ReadsWriter whenever it has been fully sequenced.
    
    The channel can be reused by calling .start() again  after a .stop(). This will however not reset the states of 
    ReadPool and ReadsWriter.
    
    Methods:
    - chan.start(t) # Start the channel at time t, channel now active
    - chan.forward(t) # Forward the channel to time t
    - chan.get_new_samples() # get new chunks of read-in-progress, concatenation of all chunks
    - chan.stop() # stop the channel, write current read until current time (last call of chan.forward(t))
    After this, the channel is clean and .start(t) can be called again
    
    A mutex protects start, stop, forward, run_mux_scan, stop_receiving, unblock.
    get_new_samples() can be called in parallel without a mutex, but modifies the stats, so they should not
    be accessed/written at the same time.

    Arguments:
        name: channel name
        read_pool: to fetch new reads as needed
        reads_writer: to write finished reads
        
    Other attributes:
        t: current simulation time, i.e. the last time to which the channel was forwarded
        t_start: time when the channel was started
        cur_elem: current element in the channel, e.g. a read; None iff the channel is not running
        stats: (cumulative) statistics of the channel until time t
        finished_elems: list of finished elements, appended to whenever save_elems is True
    
    """
    def __init__(self, name: str, read_pool: ReadPool, reads_writer: ReadsWriter, sim_params: SimParams):
        super().__init__()
        
        self.name = name
        self.read_pool = read_pool
        self.reads_writer = reads_writer
        self.sim_params = sim_params
        
        self.save_elems = False
        self.stats = None
        self.cur_elem : Union[ChannelElement, None] = None
        
        self._cur_elem_mutex = Lock()
    
    def __repr__(self):
        return f"Channel({self.name}, cur_elem={self.cur_elem}, stats={self.stats})"
    
    def start(self, t_start):
        """
        Start the channel at time t_start, starting with a ShortGap
        
        Raises:
            ChannelAlreadyRunningException: if the channel is already running
        """
        if self.is_running:
            raise ChannelAlreadyRunningException()
        
        with self._cur_elem_mutex:
            self.t_start = t_start
            self.t = t_start
            self.finished_elems = []
            self.stats = ChannelStats(n_channels=1)
            self.run_mux_scan(t_duration=0, _starting_up=True) # sets self.cur_elem
        
    def stop(self):
        """
        Stop the channel, performing necessary clean up actions (finish current element and set it to None)
        
        It will write the current read if there is one (i.e. reject it at current simulation time, i.e. last time passed to .forward()).
        Channel stats are available after calling this method.
        
        Note: The reads_writer is not flushed, so it must be done externally.
        
        Raises:
            ChannelNotRunningException: if the channel is not running
        """
        if not self.is_running:
            raise ChannelNotRunningException()
        
        with self._cur_elem_mutex:
            if isinstance(self.cur_elem, ChunkedRead):
                # reject current read
                self._write_read(self.cur_elem, end_reason=ReadEndReason.SIM_STOPPED)
            else:
                self.cur_elem.t_end = self.t
                
            self._finish_elem_in_stats(self.cur_elem)
            self.cur_elem = None
        
    @property
    def is_running(self):
        """
        Whether the channel is running (called .start(), but not .stop() yet)
        """
        return self.cur_elem is not None
    
    def is_idle(self):
        """
        If the channel is idle because no more reads are left or it is broken, but wasn't stopped yet
        """
        return isinstance(self.cur_elem, (NoReadLeftGap, ChannelBroken))
    
    def _move_to_next_elem(self, last_elem):
        """
        Helper function for .forward() to choose the next element in the channel
        Writes the current read if last_elem is a read
        
        Not thread-safe
        """
        
        t_start = last_elem.t_end
        
        # get a new read, otherwise NoReadLeftGap
        def get_new_read():
            try:
                # new read
                new_read_id, new_read_seq = self.read_pool.get_new_read(channel=self.name)
                return ChunkedRead(new_read_id, new_read_seq, t_start, t_delay=self.sim_params.gap_samplers[self.name].sample_read_start_delay(channel_stats=self.stats, random_state=self.sim_params.random_state), 
                                   read_speed=self.sim_params.bp_per_second, min_chunk_size=self.sim_params.min_chunk_size)
            except NoReadLeftException:
                # insert infinite gap
                return NoReadLeftGap(t_start)
            
        # get new gap, either short or long
        def get_new_gap():
            gap_type, gap_duration = self.sim_params.gap_samplers[self.name].sample_next_gap(channel_stats=self.stats, random_state=self.sim_params.random_state)
            if gap_type == GapSampler.GapType.Short:
                return ShortGap(t_start, t_duration=gap_duration)
            elif gap_type == GapSampler.GapType.Long:
                return LongGap(t_start, t_duration=gap_duration)
            else:
                return ChannelBroken(t_start)
            
        if isinstance(last_elem, MuxScan):
            if last_elem.elem_to_restore is None:
                new_elem = get_new_gap()
            else:
                # restore old element which is a long gap
                assert isinstance(last_elem.elem_to_restore, LongGap)
                new_elem = last_elem.elem_to_restore
                new_elem.t_start = t_start
        elif isinstance(last_elem, ChunkedRead):
            self._write_read(last_elem, end_reason=None)
            new_elem = get_new_gap()
        elif isinstance(last_elem, UnblockDelay):
            new_elem = get_new_gap()
        elif isinstance(last_elem, ShortGap):
            new_elem = get_new_read()
        elif isinstance(last_elem, LongGap):
            self.sim_params.gap_samplers[self.name].mark_long_gap_end(channel_stats=self.stats)
            new_elem = get_new_read()
        else:
            assert not isinstance(last_elem, (NoReadLeftGap, ChannelBroken)), "NoReadLeftGap has infinite length (sink state), so no next state"
            raise ValueError(f"unknown channel element type: {type(last_elem).__name__}")
        return new_elem
        
    def forward(self, t, delta=False):
        """
        Forward channel to time t
        
        Cannot go backwards in time, i.e. t >= self.t (t >= 0 if delta=True).
        
        Call .forward(t_start) at beginning to get the first read
        
        Args:
            t (float): time to forward to
            delta (bool): whether t is a delta_t (i.e. t = self.t + delta_t)
        Raises:
            ChannelNotRunningException: if channel is not running
        """
        with self._cur_elem_mutex:
            if not self.is_running:
                raise ChannelNotRunningException()
            assert self.is_running, "need to call .start(t) first"
            if delta:
                t += self.t
            assert t >= self.t, "can only forward time, not go backwards"
            
            while self.cur_elem.has_finished_by(t):
                self._update_elem_in_stats(self.cur_elem, self.t, self.cur_elem.t_end)
                self.t = self.cur_elem.t_end
                self._finish_elem_in_stats(self.cur_elem) # takes current self.t into account
                
                # important to update stats before so the gap sampling takes the updated values into account
                self.cur_elem = self._move_to_next_elem(self.cur_elem)
                self._start_elem_in_stats(self.cur_elem)
            
            self._update_elem_in_stats(self.cur_elem, self.t, t)
            self.t = t
            
            self.stats.check_consistent()
        
    ###################### functions that terminate the current element in the channel and replace it by another one #####################
    
    def run_mux_scan(self, t_duration: float, _starting_up: bool=False) -> bool:
        """
        Run a mux scan starting right now (last time forward was called), change mux scan end if one is already running.
        
        If a read is currently in the channel, it is unblocked immediately (unblock delay 0)
        If an unblock delay or short gap is happening, it is ended immediately.
        If a mux scan is already happening, it is set to end at time t_duration after the current time (self.t) of 
        the channel. This can shorten or increase the length of an existing mux scan.
        If the pore is blocked, it is split at the current location and the mux scan is inserted in 
        between keeping a reference to the second part.
        If no read is left, nothing happens.
        
        Args:
            t_duration: duration of mux scan starting from current time
            _starting_up: for internal use, used when the channel is started, mutex should already be held
        
        Returns:
            whether a read was rejected
            
        Raises:
            ChannelNotRunningException: if channel is not running
        """
        if not self.is_running and not _starting_up:
            raise ChannelNotRunningException()
        
        with self._cur_elem_mutex if not _starting_up else contextlib.nullcontext(): # starting up -> mutex already held
            assert t_duration >= 0
            elem_to_restore = None
            read_was_rejected = False
            
            if isinstance(self.cur_elem, ChunkedRead):
                # stop active read immediately
                self._write_read(self.cur_elem, end_reason=ReadEndReason.MUX_SCAN_STARTED)
                read_was_rejected = True
            elif isinstance(self.cur_elem, (UnblockDelay, ShortGap)):
                # end immediately
                self.cur_elem.t_end = self.t
            elif isinstance(self.cur_elem, LongGap):
                # split gap into two at t_split, i.e. set self to [t_start, t_split] and return a new [t_split, t_end]
                # have mux scan refer to it
                elem_to_restore = self.cur_elem.split(self.t)
            elif isinstance(self.cur_elem, MuxScan):
                # modify t_end of current mux scan, same element to restore
                self.cur_elem.t_end = self.t + t_duration
                return False # don't add MuxScan again
            elif isinstance(self.cur_elem, (NoReadLeftGap, ChannelBroken)):
                # don't do anything
                return False
            else:
                # beginning of channel, no element yet
                assert self.cur_elem is None and _starting_up, f"unknown element type {type(self.cur_elem).__name__}"
                
            # cur_elem is None when called right after start
            if self.cur_elem is not None:
                self._finish_elem_in_stats(self.cur_elem)
            
            self.cur_elem = MuxScan(self.t, t_duration=t_duration, elem_to_restore=elem_to_restore)
            self._start_elem_in_stats(self.cur_elem)
            
            return read_was_rejected
        
    def has_active_mux_scan(self) -> bool:
        return isinstance(self.cur_elem, MuxScan)
        
    def unblock(self, unblock_duration=None, end_reason=ReadEndReason.UNBLOCKED, read_id=None) -> bool:
        """
        Unblock/reject current read with given read_id.
        
        If the channel is not running or the read is not found, it is logged as a missed action.
        
        Args:
            unblock_duration: duration of unblocking delay (default: self.sim_params.default_unblock_duration)
            end_reason: reason for unblocking the read
            read_id: read id of read to unblock, if None, unblock current read
            
        Returns:
            UnblockResponse
        """
        with self._cur_elem_mutex:
            cur_elem = self.cur_elem # for thread-safety
            # add unblock delay
            if not self._write_read(self.cur_elem, end_reason=end_reason, read_id=read_id):
                # read was missed
                return UnblockResponse.MISSED
            
            self._finish_elem_in_stats(cur_elem)
            
            if unblock_duration is None:
                unblock_duration = self.sim_params.default_unblock_duration
            assert isinstance(unblock_duration, (int, float))
            
            self.cur_elem = UnblockDelay(self.t, unblock_duration, cur_elem)
            self._start_elem_in_stats(self.cur_elem) # pass in new element!
            return UnblockResponse.UNBLOCKED
    
    def _write_read(self, elem, end_reason, read_id=None) -> bool:
        """
        Finish the current read right now (possibly early) by writing it (without changing stats!)
        
        If the channel is not running, it is logged as a missed action.
        Argument order reversed to .unblock()!
        The read is not written if of length 0 (which happens due to a read delay).
        
        Returns:
            whether read was missed or not
        """
        if not isinstance(elem, ChunkedRead) or (read_id is not None and elem.full_read_id != read_id):
            # read no longer the current read
            self.stats.reads.number_rejected_missed += 1
            return False
        
        # write read up to current time t only (not necessarily full read)
        seq_record = elem.finish(self.t, end_reason=end_reason)
        seq_record.description += f" ch={self.name}"
        if DONT_WRITE_ZERO_LENGTH_READS and len(seq_record.seq) == 0:
            logger.debug(f"Read with id '{seq_record.id}' had length 0, not writing it")
        else:
            self.reads_writer.write_read(seq_record)
        
        return True
    
    
    ##################### Chunk related functions #####################
    
    def stop_receiving(self, read_id=None) -> StoppedReceivingResponse:
        """
        Stop receiving chunks from current read with read_id.
        
        If the channel is not running, it is logged as a missed action.
        
        This method should not be called in parallel from several threads (but can be called along with other methods).

        Args:
            read_id: read id of read to unblock; if None, current read

        Returns:
            True if read was stopped, False if read was not found (no longer current read)
        """
        with self._cur_elem_mutex: # for updating stats
            cur_elem = self.cur_elem # for thread-safety
            if not isinstance(cur_elem, ChunkedRead) or (read_id is not None and cur_elem.full_read_id != read_id):
                # read no longer the current read
                self.stats.reads.number_stop_receiving_missed += 1 # only method writing to this field
                return StoppedReceivingResponse.MISSED
            
            if cur_elem.stop_receiving():
                # only count if read was not already stopped
                assert self.stats.reads.cur_number == 1
                self.stats.reads.cur_number_stop_receiving += 1 # only method writing to this field
                return StoppedReceivingResponse.STOPPED_RECEIVING
            else:
                return StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING
    
    def get_new_samples(self):
        """
        Get new samples of the current read.
        
        If the read was set to stop receiving, no new samples are returned.
        
        Returns:
            Tuple of (samples, read_id, estimated_ref_len_so_far)
            If the read was set to stop_receiving, samples is ""
            If no read is active (e.g. read gap, not running), it returns ("", None, None)
        """
        
        # we are not acquiring a lock because this method will be called in parallel to "forward"
        cur_elem = self.cur_elem # for thread-safety
        if not isinstance(cur_elem, ChunkedRead):
            # also works if channel is not running
            return ("", None, None)
        chunks, read_id, estimated_ref_len_so_far = cur_elem.get_new_samples(self.t)
        self.stats.reads.number_bps_requested += len(chunks) # only method writing to this field, todo: writing racing condition
        return (chunks, read_id, estimated_ref_len_so_far)
    
    
    ##################### Channel statistics #####################
    # They all take elem as an argument of elem to make clear what they are modifying.
    # They also modify the stats, so a lock should be acquired.
    
    def _get_stats_for_elem(self, elem):
        """
        Returns object to modify in stats given current element
        """
        if isinstance(elem, ShortGap):
            return self.stats.short_gaps
        elif isinstance(elem, LongGap):
            return self.stats.long_gaps
        elif isinstance(elem, ChannelBroken):
            return self.stats.channel_broken
        elif isinstance(elem, MuxScan):
            return self.stats.mux_scans
        elif isinstance(elem, UnblockDelay):
            return self.stats.unblock_delays
        elif isinstance(elem, ChunkedRead):
            return self.stats.reads
        else:
            assert isinstance(elem, NoReadLeftGap), f"Encountered unknown element of class {elem.__class__}"
            return self.stats.no_reads_left
        
    def _start_elem_in_stats(self, elem):
        """
        Start current element in stats
        """
        self._get_stats_for_elem(elem).start()
        
    def _update_elem_in_stats(self, elem, t_from, t_to):
        """
        Add time to current element in stats
                
        Args:
            t_from, t_to: time interval [t_from, t_to] to add for current element
        """
        kwargs = {}
        if isinstance(elem, ChunkedRead):
            kwargs["nb_new_bps"] = elem.actual_seq_length(t_to) - elem.actual_seq_length(t_from) # not thread-safe
        
        self._get_stats_for_elem(elem).add_time(t_to - t_from, **kwargs)
        
    def _finish_elem_in_stats(self, elem):
        """
        Finish current element (possibly prematurely) in stats
        """
        kwargs = {}
        if isinstance(elem, ChunkedRead):
            kwargs["nb_bps_rejected"] = elem.full_seq_length() - elem.actual_seq_length(self.t)
            kwargs["stopped_receiving"] = elem.stopped_receiving
        self._get_stats_for_elem(elem).finish(**kwargs)
        
        if self.save_elems:
            self.finished_elems.append(elem)
            
    def plot(self, *args, **kwargs):
        """
        Plot channels, only plots elements recorded while save_elems was set to True
        
        Not thread-safe
        """
        return plot_channels([self], *args, **kwargs)

#%%
# time_interval: only plot elements that start or end in the time interval (if provided)
def plot_channels(channels: List[Channel], time_interval=None, ax=None, **plot_args) -> "plt.Axes":
    """
    Plot the elements in a collection of channels
    
    Requires channel elements to be recorded (i.e. channel.elem should be non-empty)
    
    Args:
        channels: list of channels to plot
        time_interval: only plot elements that start or end in the time interval (if provided)
        ax: plt axis
        plot_args: additional arguments passed when creating axes (if ax is None)
        
    Returns:
        plt axis that was plotted on
    """
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    import seaborn as sns
    
    # compute channel_positions, lines (for channel elems) and their colors
    lines = []
    colors = []
    channel_y_positions = []
    
    # create figure
    fig, ax = get_fig_ax(ax, **plot_args)

    for i, channel in enumerate(channels):
        y_pos = 0.5*(i+1)
        channel_y_positions.append(y_pos)
        for elem in channel.finished_elems + ([channel.cur_elem] if channel.cur_elem is not None else []):
            if time_interval is not None:
                if not (in_interval(elem.t_start, time_interval) or in_interval(elem.t_end, time_interval)):
                    continue
            
            if isinstance(elem, ChunkedRead):
                color = "black"
                offset = 0.03
            elif isinstance(elem, UnblockDelay):
                color = "orange"
                offset = -0.01
            elif isinstance(elem, ShortGap):
                color = "green"
                offset = 0.01
            elif isinstance(elem, LongGap):
                color = "red"
                offset = -0.03
            elif isinstance(elem, MuxScan):
                color = "purple"
                offset = -0.05
            elif isinstance(elem, NoReadLeftGap):
                color = "grey"
                offset = 0.02
            elif isinstance(elem, ChannelBroken):
                color = "blue"
                offset = 0.01
            else:
                raise TypeError(elem)
            t_end = elem.t_end
            if np.isinf(t_end):
                t_end = channel.t
                
            lines.append([(elem.t_start, y_pos + offset), (t_end, y_pos + offset)])
            colors.append(color)
        
        # element_starts = [line_segment[0][0] for line_segment in lines]
        # ax.fill_between(element_starts, y_pos-0.1, y_pos+0.1, color="grey")
    
    if len(lines) == 0:
        print("No channel elements" + ("" if time_interval is None else f" in interval {time_interval}") + "; check that save_elems was set to True")
        # return None
    
    lc = LineCollection(lines, colors=colors)
    ax.add_collection(lc)
    # add channel end times
    ax.add_collection(
        LineCollection([[(channel.t, y_pos-0.15), (channel.t, y_pos+0.15)] for (y_pos, channel) in zip(channel_y_positions, channels)], linestyle="dotted", colors="black")
    )
    ax.autoscale() # Manually adding artists doesn't rescale the plot, so we need to autoscale

    if len(lines) > 0:
        # add legend entries by creating points and making them invisible
        existing_point = lines[0][0]
        legend_elements = [
            Line2D(existing_point, existing_point, color='black', lw=2, label='read'),
            Line2D(existing_point, existing_point, color='orange', lw=2, label='unblock delay'),
            Line2D(existing_point, existing_point, color='green', lw=2, label='short gap'),
            Line2D(existing_point, existing_point, color='red', lw=2, label='long gap'),
            Line2D(existing_point, existing_point, color='purple', lw=2, label='mux scan'),
            Line2D(existing_point, existing_point, color='grey', lw=2, label='no read left'),
            Line2D(existing_point, existing_point, color='blue', lw=2, label='broken channel'),
            Line2D(existing_point, existing_point, color='black', lw=2, label='current time', linestyle="dotted"),
        ]
        ax.legend(handles=legend_elements, loc='center right')
        
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings( "ignore", module = "seaborn\..*" ) # filters MatplotlibDeprecationWarning: The legendHandles attribute was deprecated
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
            
        [elem.set_visible(False) for elem in legend_elements] # make points invisible (after creating legend)
    
    ax.set_title("Channel elements" + (" over time" if time_interval is None else f" in time interval {time_interval}"))
    ax.set_xlabel("Time (s)")

    if time_interval is not None:
        ax.set_xlim(*time_interval)
    ax.set_ylim([channel_y_positions[0]-0.5, channel_y_positions[-1]+0.5])
    # set channel names as yticks
    ax.set_yticks(channel_y_positions, [chan.name for chan in channels])
    # remove bounding box except for x-axis (bottom)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    fig.set_layout_engine(TightLayoutEngine())
    
    return ax
