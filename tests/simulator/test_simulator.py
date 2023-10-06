import itertools
import logging
from typing import Dict
import pytest
from pytest import approx
import time
import numpy as np
import matplotlib.pyplot as plt
import dill

from simreaduntil.shared_utils.timing import cur_ns_time
from simreaduntil.simulator import channel
from simreaduntil.simulator.channel import ChannelAlreadyRunningException, StoppedReceivingResponse, UnblockResponse
from simreaduntil.simulator.channel_element import ChunkedRead
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.pore_model import PoreModel
from simreaduntil.simulator.readpool import ReadPoolFromIterable
from simreaduntil.simulator.readswriter import ArrayReadsWriter
from simreaduntil.simulator.simulator import ActionType, InexistentChannelsException, ONTSimulator, ReadUntilClientFromDevice, ReadUntilDevice, assign_read_durations_to_channels, convert_action_results_to_df, plot_sim_actions, run_periodic_mux_scan_thread, stop_simulation_after_time_thread
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.simulator.utils import in_interval
from simreaduntil.usecase_helpers.utils import random_reads_gen

gen_from_list = pytest.helpers.gen_from_list
perform_random_sim_ops = pytest.helpers.perform_random_sim_ops
check_simulator_actions_agree_with_reads_writer = pytest.helpers.check_simulator_actions_agree_with_reads_writer

# cannot use autouse because it does not return its value!
@pytest.fixture()
def sim_params() -> SimParams:
    # make it fast enough so something actually happens (without making tests last too long)
    return SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.04, long_gap_length=0.05, prob_long_gap=0.02, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=100, chunk_size=20, default_unblock_duration=0.02, seed=0,
    )
    
@pytest.fixture
def simulator(sim_params) -> ONTSimulator:
    return ONTSimulator(
        read_pool=ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3), length_range=(10, 50))), 
        reads_writer=ArrayReadsWriter(),
        sim_params=sim_params
    )

def test_start_stop(simulator):
    # create simulator, start it, stop it
    # test times are in the correct range
    
    acceleration_factor = 2
    
    def test_helper(*args, **kwargs):
        start_time_max = cur_ns_time()
        simulator.start(*args, **kwargs)
        start_time_min = cur_ns_time()
        time.sleep(1.1)
        elapsed_time_min = cur_ns_time() - start_time_min
        simulator.stop()
        elapsed_time_max = cur_ns_time() - start_time_max
        assert not simulator._channels[0].is_running
        
        stats = simulator.get_channel_stats()[0]
        assert in_interval(stats.time_active + stats.no_reads_left.time_spent, ((elapsed_time_min-0.05) * acceleration_factor, elapsed_time_max * acceleration_factor)) # -0.05 since thread may not have started immediately
    
    test_helper(acceleration_factor=acceleration_factor, update_method="realtime")
    test_helper(acceleration_factor=acceleration_factor, update_method="constant")
    
    # test can be pickled
    dill.loads(dill.dumps(simulator))
    
def test_no_double_start_stop(simulator):
    # checks the simulator can't be started and stopped twice in a row, and it stays in a valid state (i.e. can still be started)
    
    assert not simulator._channels[0].is_running
    assert not simulator.stop()
    assert not simulator._channels[0].is_running
        
    simulator.start(acceleration_factor=2)
    time.sleep(0.05) # wait a bit to make sure the simulator thread has started
    assert simulator._channels[0].is_running
    
    assert not simulator.start()
    
    assert simulator._channels[0].is_running    
    time.sleep(0.2)
        
    simulator.stop()
    assert not simulator._channels[0].is_running
    
    
def test_sync_start_stop(simulator):
    simulator.sync_start()
    with pytest.raises(ChannelAlreadyRunningException):
        simulator.sync_start()
        
    simulator.get_basecalled_read_chunks()
    for (i, delta_t) in enumerate([0.2, 0.1, 0.5, 1, 2, 10, 0.05, 0.4]):
        simulator.sync_forward(delta_t, delta=True)
        if i % 3 == 2:
            simulator.get_basecalled_read_chunks()
    simulator.sync_stop()

def test_get_basecalled_chunks():
    # test channel subset for chunks, batch size
    
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=0.5, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=0.2, seed=0,
    )
    
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGC"), ("read2", "TTTTAC"), ("read3", "TTTTAAAACCCAAACTTTACCA"), ("read4", "TCTTAAAACCTTA"))))
    simulator = ONTSimulator(
        read_pool=read_pool,
        reads_writer=ArrayReadsWriter(),
        sim_params = sim_params
    )
    simulator.save_elems = True
    eps = 1e-5
    
    simulator.sync_start(0)
    simulator.sync_forward(0.9)
    chunks = list(simulator.get_basecalled_read_chunks())
    assert sorted(chunks) == sorted([(1, 'read1', 'AAAA', 'noquality', 4), (2, 'read2', 'TTTT', 'noquality', 4)])
    
    chunks = list(simulator.get_basecalled_read_chunks())
    assert len(chunks) == 0
    
    simulator.sync_forward(1.2+eps)
    chunks = list(simulator.get_basecalled_read_chunks())
    assert chunks == [(1, 'read1', 'GGGG', 'noquality', 8)]
    
    simulator.sync_forward(1.4+eps) # to force channel 1 to get read3
    simulator.sync_forward(2.2+eps)
    with pytest.raises(InexistentChannelsException):
        # channels are 1-based
        list(simulator.get_basecalled_read_chunks(channel_subset=[0]))
    chunks = list(simulator.get_basecalled_read_chunks(channel_subset=[1]))
    assert chunks == [(1, 'read4', 'TCTT', 'noquality', 4)]
    chunks = list(simulator.get_basecalled_read_chunks(channel_subset=[2]))
    assert chunks == [(2, 'read3', 'TTTTAAAA', 'noquality', 8)]
    
    simulator.sync_forward(2.7+eps)
    chunks = list(simulator.get_basecalled_read_chunks(batch_size=1))
    assert len(chunks) == 1
    
    # simulator.plot_channels(); import matplotlib.pyplot as plt; plt.show()
    
    simulator.sync_stop()
    
    chunks = list(simulator.get_basecalled_read_chunks(batch_size=1))
    assert len(chunks) == 0
    
def test_get_raw_chunks(shared_datadir):
    seq = "ACCCTTTGGGCCGG"
    k = 6
    signals_per_bp = 2
    pore_filename = shared_datadir / "dummy_pore_model.csv"
    eps = 1e-5
    
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=0.5, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=0.2, seed=0, pore_model=PoreModel(pore_filename, signals_per_bp=signals_per_bp)
    )
    
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", seq), ("read2", "TTTTAC"))))
    simulator = ONTSimulator(
        read_pool=read_pool,
        reads_writer=ArrayReadsWriter(),
        sim_params = sim_params
    )
    
    simulator.sync_start(0)
    
    simulator.sync_forward(0.4+0.9+eps)
    chunks = list(simulator.get_raw_chunks(channel_subset=[1])) # get 2 chunks for channel 1
    raw_signal = chunks[0][2]
    assert len(raw_signal) == (2*4 - k + 1) * signals_per_bp
    
def test_synchronous_sim_is_deterministic(sim_params, channel_write_zero_length_reads):
    # test that the synchronous simulator produces the same results when run twice, only for "constant" update method
    # async mode is not deterministic because it is an ongoing loop (runs in the background at unregular intervals)
    
    async_mode = False 
    
    def run_sim(seed):
        reads_writer = ArrayReadsWriter()
        
        simulator = ONTSimulator(
            read_pool=ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3), length_range=(10, 50))), 
            reads_writer=reads_writer,
            sim_params=sim_params
        )
        simulator.save_elems = True
        
        acceleration_factor = 10
        simulator.start(acceleration_factor=acceleration_factor, update_method="constant", seed=seed) if async_mode else simulator.sync_start(0)
        
        nb_actions = perform_random_sim_ops(simulator, random_state=np.random.default_rng(seed), acceleration_factor=acceleration_factor, nb_iterations=100)
        
        simulator.stop() if async_mode else simulator.sync_stop()
        
        check_simulator_actions_agree_with_reads_writer(simulator, nb_actions)
        
        return (reads_writer.reads, nb_actions, simulator)
    
    for seed in range(2):
        # a1 = run_sim(seed=seed)
        # a2 = run_sim(seed=seed)
        assert run_sim(seed=seed)[:2] == run_sim(seed=seed)[:2]


@pytest.mark.parametrize("async_mode", [(False,), (True,)])
def test_random_ops_synchronous(simulator, async_mode, channel_write_zero_length_reads):
    # perform random operations and check that the written reads correspond to the actions that were performed
    
    # make all elements roughly the same length, reads have length 8-17 -> take about 1.2s
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=1.2, long_gap_length=1.2, prob_long_gap=0.35, time_until_blocked=200, read_delay=0) for i in range(2)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=1.2, seed=0,
    )
    
    # apply random operations, check that reads are correct
    all_read_ids = []
    def get_read_and_save_id():
        for (read_id, read_seq) in random_reads_gen(random_state=np.random.default_rng(2)):
            all_read_ids.append(read_id)
            yield read_id, read_seq
           
    simulator = ONTSimulator(
        read_pool=ReadPoolFromIterable(get_read_and_save_id()), 
        reads_writer=ArrayReadsWriter(),
        sim_params=sim_params
    )
    simulator.save_elems = True
    
    cum_nb_actions = None
    for j in range(3):
        print(f"j={j}", end=", ")
        
        t_start = 2
        simulator.start(acceleration_factor=20) if async_mode else simulator.sync_start(t_start)
        
        nb_actions = perform_random_sim_ops(simulator, random_state=np.random.default_rng(j), async_mode=async_mode, acceleration_factor=20, nb_iterations=200)
        
        if cum_nb_actions is None:
            cum_nb_actions = nb_actions
        else:
            cum_nb_actions = {k: cum_nb_actions[k] + nb_actions[k] for k in cum_nb_actions}
        
        if nb_actions["forward"] > 0:
            pass
            # print(1)
            # ax = simulator.plot_channels(); import matplotlib.pyplot as plt; plt.show(block=True); save_fig_and_pickle(ax.figure, "simulator_example.png");  #note: for README.md figure, zoom into appropriate region, then save the figure, keep it commented out
            
        simulator.stop() if async_mode else simulator.sync_stop()
        check_simulator_actions_agree_with_reads_writer(simulator, cum_nb_actions, check_nb_rejected=(j==0))
        
        stats = simulator.get_channel_stats()[0]
        # there can be a difference due to the last iteration until the simulator is fully stopped
        assert stats.time_active + stats.no_reads_left.time_spent + stats.channel_broken.time_spent == approx(simulator._channels[0].t - (0 if async_mode else t_start))    
    
def test_realtime(channel_write_zero_length_reads):
    # test that the simulator can run in real time (i.e. end time is at least 0.95 * real time)
    # can be used to determine the optimal acceleration factor that does not cause too much delay
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=1.3, prob_long_gap=0.1, time_until_blocked=np.inf, read_delay=0) for i in range(512)},
        bp_per_second=450, chunk_size=200, default_unblock_duration=1.4, seed=0,
    )
    
    # pre-generate reads to make sure that this is not responsible for the delay
    # reads_gen = random_reads_gen(random_state=np.random.default_rng(3), length_range=(500, 5000))
    
    simulator = ONTSimulator(
        # read_pool=ReadPoolFromIterable(gen_from_list(itertools.islice(reads_gen, 10000))),
        # read_pool=ReadPoolFromIterable(reads_gen), 
        read_pool=ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3), length_range=(500, 5000))),
        reads_writer=ArrayReadsWriter(),
        sim_params=sim_params
    )
    
    acceleration_factor = 5 # depends on computer load
    simulator.start(acceleration_factor=acceleration_factor)
    t_start = cur_ns_time()
    nb_actions = perform_random_sim_ops(simulator, random_state=np.random.default_rng(4), acceleration_factor=acceleration_factor, nb_iterations=500)
    simulator.stop() # if it fails, check enough reads were pre-generated above
    check_simulator_actions_agree_with_reads_writer(simulator, nb_actions)
    
    time_elapsed = cur_ns_time() - t_start
    
    stats = simulator.get_channel_stats()[0]
    simulated_time = stats.time_active + stats.no_reads_left.time_spent
    min_sim_time_required = (time_elapsed-0.05) * acceleration_factor * 0.95 # 0.05 for sim thread startup time
    print(f"Simulated time: {simulated_time}, min sim time required: {min_sim_time_required}, real time: {time_elapsed}")
    assert simulated_time >= min_sim_time_required, f"cannot really keep up with real time at acceleration factor {acceleration_factor}"

def test_stop_simulation_after_time_thread(simulator):
    simulator.start()
    stop_simulation_after_time_thread(simulator, 1.0).start()
    time.sleep(1.5) # 1.15 sometimes fails on CI
    assert not simulator.is_running
    
def test_run_periodic_mux_scan_thread(simulator):
    simulator.start()
    mux_scan_thread = run_periodic_mux_scan_thread(simulator, period=1.0, scan_duration=0.1, acceleration_factor=3)
    mux_scan_thread.start()
    time.sleep(4.3/3)
    simulator.stop()
    
    # mux scan at 0 and [0.0, 1.0, 2.0, 3.0, 4.0] roughly
    assert simulator.get_channel_stats()[0].mux_scans.finished_number == 5
    assert simulator.get_channel_stats()[1].mux_scans.finished_number == 5
    
    time.sleep(1.0) # check that the thread does not crash because the simulation is no longer running
    mux_scan_thread.raise_if_error()
    
    assert simulator.get_channel_stats()[0].mux_scans.finished_number == 5
    assert simulator.get_channel_stats()[1].mux_scans.finished_number == 5
    
def test_readuntil_fromdevice():
    # tests the readuntil client
    
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=0.5, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=0.2, seed=0,
    )
    
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGC"), ("read2", "TTTTACCTTACC"), ("read3", "TTTTAAAACCCAAACTTTACCA"), ("read4", "TCTTAAAACCTTA"))))
    simulator = ONTSimulator(
        read_pool=read_pool,
        reads_writer=ArrayReadsWriter(),
        sim_params=sim_params
    )
    simulator.save_elems = True
    eps = 1e-5
    
    ru_client = ReadUntilClientFromDevice(simulator)
    
    simulator.sync_start(-0.4)
    simulator.sync_forward(0.5)
    
    chunks = list(ru_client.get_basecalled_read_chunks())
    assert len(chunks) == 2
    
    with pytest.raises(InexistentChannelsException):
        ru_client.stop_receiving_batch([(0, "read1")])
    responses = ru_client.stop_receiving_batch([(1, "read1"), (2, "inexistent")])
    assert responses == [StoppedReceivingResponse.STOPPED_RECEIVING, StoppedReceivingResponse.MISSED]
    assert ru_client.stop_receiving_batch([(1, "read1")]) == [StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING]
    
    simulator.sync_forward(0.8+eps)
    
    chunks = list(ru_client.get_basecalled_read_chunks())
    assert len(chunks) == 1
    
    assert ru_client.stop_receiving_batch([(2, "read2")]) == [StoppedReceivingResponse.STOPPED_RECEIVING]
    
    assert ru_client.unblock_read_batch([(1, "read1"), (2, "read2"), (1, "inexistent")]) == [True, True, False]
    
    simulator.sync_stop()
    
def test_get_action_results():
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=0.5, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=0.2, seed=0,
    )
    
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGC"), ("read2", "TTTTAC"), ("read3", "TTTTAAAACCCAAACTTTACCA"), ("read4", "TCTTAAAACCTTA"))))
    simulator = ONTSimulator(
        read_pool=read_pool,
        reads_writer=ArrayReadsWriter(),
        sim_params = sim_params
    )
    simulator.save_elems = True
    
    simulator.sync_start(-0.4)
    simulator.sync_forward(0.5)
    simulator.stop_receiving_read(1, "read1")
    simulator.stop_receiving_read(2, "inexistent")
    exp_action_results = [
        ("read1", 0.5, 1, ActionType.StopReceiving, StoppedReceivingResponse.STOPPED_RECEIVING),
        ("inexistent", 0.5, 2, ActionType.StopReceiving, StoppedReceivingResponse.MISSED)
    ]
    assert simulator.get_action_results(clear=False) == exp_action_results
    plot_sim_actions(convert_action_results_to_df(exp_action_results), close_figures=True)
    
    assert simulator.get_action_results() == exp_action_results
    assert len(simulator.get_action_results()) == 0
    simulator.sync_forward(1.1+1e-4) # needed so that channel 2 gets read3
    simulator.sync_forward(1.3+1e-4)
    simulator.unblock_read(2, "read3")
    simulator.unblock_read(2, "inexistent")
    assert simulator.get_action_results() == [
        ("read3", 1.3+1e-4, 2, ActionType.Unblock, UnblockResponse.UNBLOCKED),
        ("inexistent", 1.3+1e-4, 2, ActionType.Unblock, UnblockResponse.MISSED)
    ]
    assert simulator._channels[0].cur_elem.full_read_id == "read4"
    
    simulator.sync_stop() # does not perform a user action
    assert len(simulator.get_action_results()) == 0

def test_assign_read_durations_to_channels():
    # too few read durations
    assert assign_read_durations_to_channels(
        [np.array([0.1, 0.2, 0.3]), np.array([0.1, 0.2, 0.3])+1], 
        ["ch1", "ch2", "ch3"]
    ) == {"ch1": approx(np.array([0.1, 0.2, 0.3])), "ch2": approx(np.array([0.1, 0.2, 0.3])+1), "ch3": approx(np.array([0.1, 0.2, 0.3]))}
    
    # too many read durations
    assert assign_read_durations_to_channels(
        [np.array([0.1, 0.2, 0.3]), np.array([0.1, 0.2, 0.3])+1, np.array([0.1, 0.2, 0.3])+2, np.array([0.1, 0.2, 0.3])+3], 
        ["ch1", "ch2", "ch3"]
    ) == {"ch1": approx(np.array([0.1, 0.2, 0.3])), "ch2": approx(np.array([0.1, 0.2, 0.3])+1), "ch3": approx(np.array([0.1, 0.2, 0.3])+2)}