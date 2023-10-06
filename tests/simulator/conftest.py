# helper functions for tests
# for some reason, pytest only makes fixtures defined in conftest available, not regular functions
# use pytest-helpers-namespace
# disadvantage: pylint does not recognize it
# alternative:
# use "Factories as fixtures", autouse does not work (error: fixtures are meant to be called directly), so we have to list it as parameters explicitly and write
# things like: gen_from_list = pytest.fixture()(_gen_from_list)

from contextlib import contextmanager
import time
from typing import Dict
import numpy as np
import pytest
from simreaduntil.shared_utils.dna import get_random_DNA_seq
from simreaduntil.simulator import channel
from simreaduntil.simulator.channel import StoppedReceivingResponse, UnblockResponse
from simreaduntil.simulator.channel_element import ReadDescriptionParser
from simreaduntil.simulator.readswriter import ArrayReadsWriter
from simreaduntil.simulator.simulator import ActionType, ONTSimulator, ReadUntilDevice

@pytest.fixture()
def channel_write_zero_length_reads():
    # temporarily write reads of length zero
    # use fixture in case the test fails
    prev_val = channel.DONT_WRITE_ZERO_LENGTH_READS
    channel.DONT_WRITE_ZERO_LENGTH_READS = False # so that we see all reads passed to the readswriter
    yield None
    channel.DONT_WRITE_ZERO_LENGTH_READS = prev_val

@pytest.helpers.register
def gen_from_list(lst):
    """
    Create a generator from a list
    """
    return (x for x in lst)

@pytest.helpers.register
def dict_different_keys(d_expected, d2):
    """
    For two dicts with approximately equal keys, return list of keys where values are different
    """
    assert d_expected.keys() == d2.keys(), f"keys not equal: {set(d_expected.keys()).symmetric_difference(d2.keys())}"
    
    from pytest import approx
    unequal_keys = [k for k in d_expected.keys() if d_expected[k] != approx(d2[k])]
    print(f"unequal keys: {unequal_keys}")
    print(f"expected vals: {[d_expected[k] for k in unequal_keys]}")
    print(f"  actual vals: {[d2[k] for k in unequal_keys]}")
    return unequal_keys


############################################
# helper functions for the simulator
############################################

@contextmanager
def WriteZeroLengthReadsSetting():
    # Context manager to temporarily set DONT_WRITE_ZERO_LENGTH_READS
    # so that we see all reads passed to the readswriter
    prev_val = channel.DONT_WRITE_ZERO_LENGTH_READS
    channel.DONT_WRITE_ZERO_LENGTH_READS = False
    yield None
    channel.DONT_WRITE_ZERO_LENGTH_READS = prev_val

@pytest.helpers.register
def perform_random_sim_ops(simulator: ReadUntilDevice, random_state: np.random.Generator, async_mode=True, acceleration_factor=1, nb_iterations=100) -> Dict[str, int]:
    """
    Performs random operations on the simulator
    
    It either performs a readuntil action (unblock, stop_receiving, get_chunks), a mux scan or a forward (just sleeping in synchronous mode).
    
    Use the context manager WriteZeroLengthReadsSetting() or channel_write_zero_length_reads fixture to write zero length reads.
    This cannot be integrated into this function since this function does not stop the simulator (and a zero length read may occur when the simulator is stopped).
    
    Args:
        simulator: the simulator to perform operations on
        random_state: the state to use to sample operations to perform
        async_mode: whether to run in async or synchronous mode
        acceleration_factor: acceleration factor at which the simulation is running
        nb_iterations: number of iterations to perform
    
    Returns: number of actions of each type (dict)
    """
    last_read_id_per_channel = {chan: None for chan in range(1, 1 + simulator.n_channels)} # last read id per channel received from get_basecalled_read_chunks
    def get_random_chan_and_read_id():
        chan = random_state.choice(simulator.n_channels) + 1
        # read_id = None means current read (if any)
        read_id = random_state.choice([last_read_id_per_channel[chan], "inexistent"], p=[0.9, 0.1])
        
        return {"read_channel": chan, "read_id": read_id}
    
    forward_sim = lambda delta_t: time.sleep(delta_t / acceleration_factor) if async_mode else simulator.sync_forward(delta_t * acceleration_factor, delta=True)
    
    
    nb_actions = {"forward": 0, "stop_receiving": 0, "user_unblock": 0, "mux_scan_nb_unblocked": 0, "get_new_chunks": 0, "run_mux_scan": 0}
    i = 0
    last_printed_time = None
    while i < nb_iterations:
        i += 1
        if last_printed_time is None or time.time() - last_printed_time > 1:
            last_printed_time = time.time()
            print(f"iteration {i}", end="\r")
        
        weights = np.array([0.9, 0.044, 0.044, 0.01])
        u = random_state.choice(["forward", "stop_receiving", "user_unblock", "run_mux_scan"], p=weights/weights.sum())
        nb_actions[u] += 1
        if u == "forward":
            forward_sim(random_state.uniform(0.3, 0.6))
        elif u == "stop_receiving":
            simulator.stop_receiving_read(**get_random_chan_and_read_id())
        elif u == "user_unblock":
            simulator.unblock_read(**get_random_chan_and_read_id(), unblock_duration=random_state.uniform(0.3, 0.6))
        elif u == "run_mux_scan":
            t_duration = random_state.uniform(0.3, 0.6)
            nb_actions["mux_scan_nb_unblocked"] += simulator.run_mux_scan(t_duration=t_duration)
            # not exactly the same as sum(isinstance(chan.cur_elem, ChunkedRead) for chan in simulator._channels) in async mode because channels can be forwarded while mux scan is started 
            
            forward_sim(t_duration + 0.01)
        
        if random_state.uniform() < 0.2:
            for (chan, read_id, *_) in simulator.get_basecalled_read_chunks():
                last_read_id_per_channel[chan] = read_id
    
    return nb_actions

@pytest.helpers.register
def check_simulator_actions_agree_with_reads_writer(simulator: ONTSimulator, nb_actions, check_nb_rejected=True):
    """
    Check that the action results performed on the simulator agree with the reads writer
    
    Args:
        nb_actions: cumulative number of actions performed on the simulator, not necessarily all performed successfully, as an upper bound
        check_nb_rejected: check number of rejected basepairs from stats matches the one in the reads writer, only true until before the second start (since stats is reset, readswriter usually not)
    """
    action_results = simulator.get_action_results(clear=False)
    reads_writer = simulator._reads_writer
    assert isinstance(reads_writer, ArrayReadsWriter)
    written_stopped_receiving = sum(["stopped_receiving" in read_info[2] for read_info in reads_writer.reads])
    written_user_unblocked = sum(["user_unblocked" in read_info[2] for read_info in reads_writer.reads]) # only count user unblocks
    written_mux_unblocked = sum(["mux_scan_unblocked" in read_info[2] for read_info in reads_writer.reads]) # mux unblocked
    
    assert written_stopped_receiving <= nb_actions["stop_receiving"]
    assert written_user_unblocked <= nb_actions["user_unblock"]
    assert written_mux_unblocked <= nb_actions["mux_scan_nb_unblocked"]
    
    # only count the first stopped_receiving
    # (id, time, channel, action_type, action_result)
    assert sum(res[3] == ActionType.StopReceiving and res[4] == StoppedReceivingResponse.STOPPED_RECEIVING for res in action_results) == written_stopped_receiving
    sum(res[3] == ActionType.Unblock and res[4] == UnblockResponse.UNBLOCKED for res in action_results) == written_user_unblocked
    
    if check_nb_rejected:
        # check number of rejected basepairs is correct by parsing reads_writer description and comparing to simulator combined stats
        assert not simulator.is_running # otherwise race condition with forward
        written_number_bps_rejected = sum(ReadDescriptionParser(desc).full_seqlen - len(seq) for (_, seq, desc) in reads_writer.reads)
        assert simulator.get_channel_stats(combined=True).reads.number_bps_rejected == written_number_bps_rejected, "check that you have used the fixture channel_write_zero_length_reads"
