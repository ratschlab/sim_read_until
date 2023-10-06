import pytest
from pytest import approx
import numpy as np
import matplotlib.pyplot as plt
import dill

from simreaduntil.simulator import channel
from simreaduntil.simulator.channel import Channel, ChannelAlreadyRunningException, ChannelNotRunningException, StoppedReceivingResponse, plot_channels
from simreaduntil.simulator.channel_element import LongGap, ChunkedRead
from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.readpool import ReadPoolFromIterable
from simreaduntil.simulator.readswriter import ArrayReadsWriter, SingleFileReadsWriter
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.usecase_helpers.utils import random_reads_gen

gen_from_list = pytest.helpers.gen_from_list

# cannot use autouse because it does not return its value!
@pytest.fixture()
def sim_params() -> SimParams:
    return SimParams(
        gap_samplers={"channel1": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=1.4, seed=0,
    )
    
@pytest.fixture()
def channel_stats() -> ChannelStats:
    # returns channel stats with 1 channel and initial mux scan and a short gap
    channel_stats = ChannelStats(n_channels=1)
    channel_stats.mux_scans.start()
    channel_stats.mux_scans.finish()
    channel_stats.short_gaps.start_and_add_time(0.4)
    channel_stats.short_gaps.finish()
    return channel_stats


eps = 1e-8 # small delay (to avoid issues with rounding errors when geting chunks up to time <= t)
    
def test_channel_readwriting(sim_params):

    # check start and stop, stop should write to file
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    assert not chan.is_running
    chan.start(1-0.4) # first mux scan followed by short_gap, so shifting by 0.4
    assert chan.is_running
    chan.forward(1+ 1.+eps)
    assert reads_writer.reads == []
    assert chan.is_running
    chan.stop()
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGGGGCC")]
    assert not chan.is_running
    
    # check channel can be started again, at smaller time
    chan.start(-2-0.4)
    assert chan.is_running
    chan.forward(-2+0.6+eps)
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGGGGCC")]
    assert chan.is_running
    chan.stop()
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGGGGCC"), ("read2", "TTTTAA")]
    
    # test can be pickled
    dill.loads(dill.dumps(chan))
    
def test_channel_readwriting_zero_length(sim_params, channel_write_zero_length_reads):
    # this test may fail when running tests in parallel (if same python module is used across tests)
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(1.0)
    chan.forward(1.0 + 0.4 + eps)
    chan.stop()
    
    assert [x[:2] for x in reads_writer.reads] == [("read1", "")]
    
# check that channel stats are correctly updated
def test_channel_stats(sim_params, channel_stats: ChannelStats):
    # reads of length 14, 10, 7
    t_start = 2
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"), ("read3", "AATTTCT"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(t_start)
    
    chan.forward(t_start + 1.1 + eps)
    channel_stats.reads.start_and_add_time(0.7 + eps, nb_new_bps=7)
    assert chan.stats == channel_stats
    
    chan.get_new_chunks()
    channel_stats.reads.number_bps_requested = 4
    assert chan.stats == channel_stats
    
    chan.stop_receiving()
    chan.get_new_chunks() # no chunks
    channel_stats.reads.cur_number_stop_receiving = 1
    assert chan.stats == channel_stats
    
    chan.stop_receiving() # again should have no effect
    channel_stats.reads.cur_number_stop_receiving = 1
    assert chan.stats == channel_stats
    
    chan.forward(t_start + 1.7 + eps)
    channel_stats.reads.add_time(0.6, nb_new_bps=6)
    assert chan.stats == channel_stats
    
    chan.forward(t_start + 1.9 + eps)
    channel_stats.reads.add_time_and_finish(0.1, nb_new_bps=1, stopped_receiving=True)
    channel_stats.short_gaps.start_and_add_time(0.1)
    assert chan.stats == channel_stats
    
    chan.stop_receiving("inexistent")
    channel_stats.reads.number_stop_receiving_missed = 1
    assert chan.stats == channel_stats
    
    chan.forward(t_start + 2.3 + eps)
    channel_stats.short_gaps.add_time_and_finish(0.3)
    channel_stats.reads.start_and_add_time(0.1, nb_new_bps=1)
    assert chan.stats == channel_stats
    
    # unblock
    chan.unblock(0.3)
    chan.forward(t_start + 2.3 + 0.15)
    channel_stats.reads.add_time_and_finish(0, nb_new_bps=0, stopped_receiving=False, nb_bps_rejected=9)
    channel_stats.unblock_delays.start_and_add_time(0.15)
    assert chan.stats == channel_stats
    
    # move until all reads are depleted
    chan.forward(t_start + 10.3 + eps)
    channel_stats.unblock_delays.add_time_and_finish(0.15)
    channel_stats.short_gaps.add_full(0.4)
    channel_stats.reads.add_full(0.7, nb_new_bps=7, stopped_receiving=False)
    channel_stats.short_gaps.add_full(0.4)
    channel_stats.no_reads_left.start_and_add_time(10.3 + eps - channel_stats.time_active)
    assert chan.stats == channel_stats
    
    chan.stop()
    channel_stats.no_reads_left.finish()
    assert channel_stats.n_channels_running == 0
    assert chan.stats == channel_stats
    
def test_unblocking(sim_params, channel_stats):
    t_start = 2
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.save_elems = True
    chan.start(t_start)
    chan.forward(t_start + 1.0 + eps)
    
    channel_stats.reads.start_and_add_time(0.6 + eps, nb_new_bps=6)
    assert chan.stats == channel_stats
    
    assert chan.unblock(0.3)
    channel_stats.reads.finish(stopped_receiving=False, nb_bps_rejected=14-6)
    channel_stats.unblock_delays.start()
    assert chan.stats == channel_stats
    
    chan.forward(t_start + 1.2 + eps)
    channel_stats.unblock_delays.add_time(0.2)
    assert chan.stats == channel_stats
    
    assert not chan.unblock(0.3), "no read to unblock"
    
# unblock on read that was set to stop_receiving
def test_stopreceiving_then_unblock(sim_params, channel_stats: ChannelStats):
    t_start = 2
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(t_start)
    chan.forward(t_start + 1.0 + eps)
    chan.get_new_chunks()
    assert chan.stop_receiving() == StoppedReceivingResponse.STOPPED_RECEIVING
    assert chan.stop_receiving() == StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING
    
    channel_stats.reads.start_and_add_time(0.6 + eps, nb_new_bps=6)
    channel_stats.reads.cur_number_stop_receiving += 1
    channel_stats.reads.number_bps_requested += 4
    assert chan.stats == channel_stats
    
    assert chan.unblock(0.3)
    
    channel_stats.reads.finish(stopped_receiving=True, nb_bps_rejected=14-6)
    channel_stats.unblock_delays.start()
    assert chan.stats == channel_stats
    
    assert not chan.unblock(0.3), "no read to unblock"
    assert chan.stop_receiving() == StoppedReceivingResponse.MISSED, "no read to unblock"
    channel_stats.reads.number_rejected_missed += 1
    channel_stats.reads.number_stop_receiving_missed += 1
    assert chan.stats == channel_stats
    
    chan.forward(t_start + 1.2 + eps)
    channel_stats.unblock_delays.add_time(0.2)
    assert chan.stats == channel_stats
    
    chan.stop()
    
def test_channel_restart(sim_params, channel_stats: ChannelStats):
    t_start = 2
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(t_start)
    with pytest.raises(ChannelAlreadyRunningException):
        chan.start(t_start + 1)
        
    chan.forward(t_start + 1.0 + eps)
    channel_stats.reads.start_and_add_time(0.6 + eps, nb_new_bps=6)
    assert chan.stats == channel_stats
    
    # check restart by not replacing Channel, only read_pool and reads_writer
    chan.stop()
    
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan.read_pool = read_pool
    chan.reads_writer = reads_writer
    
    chan.start(t_start)
    chan.forward(t_start + 1.0 + eps)
    assert chan.stats == channel_stats
    
    chan.stop()
    with pytest.raises(ChannelNotRunningException):
        chan.stop()
    
    channel_stats.reads.finish(stopped_receiving=False, nb_bps_rejected=14-6)
    assert channel_stats.n_channels_running == 0
    assert chan.stats == channel_stats
    
def test_get_new_chunks(sim_params):
    # reads of length 14, 10, 10, 10
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"), ("read3", "TTTTAAAACC"), ("read4", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    # chan.save_elems = True # for plotting afterwards
    chan.start(-0.4)
    chan.forward(0.3 + eps)
    assert [x[:2] for x in reads_writer.reads] == []
    assert chan.get_new_chunks()[:2] == ("", "read1")
    chan.forward(0.4 + eps)
    assert chan.get_new_chunks()[:2] == ("AAAA", "read1")
    assert chan.get_new_chunks()[:2] == ("", "read1") # already returned chunks
    assert [x[:2] for x in reads_writer.reads] == []
    chan.forward(1.3 + eps)
    assert chan.get_new_chunks()[:2] == ("GGGGCCCC", "read1")
    assert [x[:2] for x in reads_writer.reads] == []
    chan.forward(1.4 + eps)
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGGGGCCCCTT")]
    
    # read entire second read
    chan.forward(0.4 + 2.4 + eps) # 14+10
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGGGGCCCCTT"), ('read2', 'TTTTAAAACC')]

    assert chan.get_new_chunks()[:2] == ("", None)
    chan.forward(2*0.4 + 3.1 + eps) # 14+10+7
    assert chan.get_new_chunks()[:2] == ("TTTT", "read3")
    
    chan.stop()
    chan.get_new_chunks() # no exception
    
    # ax = plot_channels([chan], time_interval=[0.9, 6.2], figsize=(6, 2))
    # ax.figure.show()

def test_read_stop_receiving(sim_params):
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(2-0.4)
    chan.forward(2.6 + eps)
    assert chan.get_new_chunks()[:2] == ("AAAA", "read1")
    assert chan.stop_receiving() == StoppedReceivingResponse.STOPPED_RECEIVING
    assert chan.stop_receiving() == StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING
    assert chan.get_new_chunks()[:2] == ("", "read1")
    
    assert chan.stop_receiving(read_id="inexistent") == StoppedReceivingResponse.MISSED
    
    chan.forward(3.5 + eps)
    assert chan.stop_receiving() == StoppedReceivingResponse.MISSED
    assert chan.get_new_chunks()[:2] == ("", None), "in a gap"

def test_mux_scan(sim_params, channel_stats):
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    
    chan.start(2-0.4)
    chan.forward(2.6 + eps)
    
    channel_stats.reads.start_and_add_time(0.6 + eps, nb_new_bps=6)
    assert chan.stats == channel_stats
    
    assert not chan.has_active_mux_scan()
    chan.run_mux_scan(t_duration=1.0)
    assert chan.has_active_mux_scan()
    
    channel_stats.reads.finish(stopped_receiving=False, nb_bps_rejected=14-6)
    channel_stats.mux_scans.start()
    assert chan.stats == channel_stats
    
    chan.run_mux_scan(t_duration=1.2) # modify end of mux scan
    assert chan.stats == channel_stats
    
    chan.forward(3.1)
    channel_stats.mux_scans.add_time(0.5-eps)
    assert chan.stats == channel_stats
    
    chan.run_mux_scan(t_duration=0.3) # shorten again
    assert chan.stats == channel_stats
    
    # move past mux scan
    chan.forward(3.9 + 2*eps)
    channel_stats.mux_scans.add_time_and_finish(0.3+eps)
    channel_stats.short_gaps.add_full(0.4)
    channel_stats.reads.start_and_add_time(0.1+eps, nb_new_bps=1)
    assert chan.stats == channel_stats
    
    chan.stop()
    
    assert [x[:2] for x in reads_writer.reads] == [("read1", "AAAAGG"), ("read2", "T")]
    
    with pytest.raises(ChannelNotRunningException):
        chan.run_mux_scan(t_duration=1.0)
    
def test_poreblockage_continues_after_mux_scan():
    # long pore blockage: check that it continues after the mux scan
    
    sim_params = SimParams(
        gap_samplers={"channel1": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.4, prob_long_gap=0.9, time_until_blocked=np.inf, read_delay=0)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=1.4, seed=0,
    )
    read_pool = ReadPoolFromIterable(gen_from_list((("read1", "AAAAGGGGCCCCTT"), ("read2", "TTTTAAAACC"))))
    reads_writer = ArrayReadsWriter()
    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    chan.save_elems = True
    
    chan.start(10)
    while not isinstance(chan.cur_elem, LongGap):
        chan.forward(0.3, delta=True)
    assert chan.cur_elem.t_duration == 10.4
    blocked_start = chan.cur_elem.t_start
    channel_time = chan.t
    chan.forward(2, delta=True)
    chan.run_mux_scan(t_duration=1.0)
    chan.run_mux_scan(t_duration=1.5)
    chan.forward(3, delta=True)
    chan.run_mux_scan(t_duration=1.2)
    chan.run_mux_scan(t_duration=0.8)
    chan.forward(1.6, delta=True)
    
    assert isinstance(chan.cur_elem, LongGap)
    assert chan.cur_elem.t_start == channel_time + 2 + 3 + 0.8
    assert chan.cur_elem.t_end == blocked_start + 10.4 + 1.5 + 0.8
    
    # from matplotlib import pyplot as plt; chan.plot(); plt.show()
    
def test_plotting(sim_params):
    
    sim_params.set(gap_samplers={f"channel{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0.5, time_until_blocked=np.inf, read_delay=0) for i in range(1, 4)})
    
    read_pool = ReadPoolFromIterable(random_reads_gen())
    reads_writer = ArrayReadsWriter()
    random_state = np.random.default_rng(2)
    
    # plot several channels
    channels = []
    for i in range(1, 4):
        chan = Channel(f"channel{i}", read_pool, reads_writer, sim_params=sim_params)
        chan.save_elems = True # for plotting afterwards
        t_start = 5
        chan.start(t_start)
        for delta_t in random_state.uniform(0.1, 20, size=(50,)):
            chan.forward(delta_t, delta=True)
            if isinstance(chan.cur_elem, LongGap):
                # print(chan.cur_elem)
                pass
            else:
                chunks = chan.get_new_chunks()[0]
                # if len(chunks) > 0:
                #     print(f"{delta_t}: {chunks}")
        channels.append(chan)
    
    ax = chan.plot(); plt.close(ax.figure)
    chan.stop()
        
    ax = plot_channels(channels, time_interval=[50, 350], figsize=(6, 2)); plt.close(ax.figure)
    
    # ax.figure.show()
    # save_fig_and_pickle(ax.figure, "channel_occupation.pdf")

def test_channel_normal_operation(sim_params):
    reads_writer = SingleFileReadsWriter()
    read_pool = ReadPoolFromIterable(random_reads_gen())

    chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
    t_start = 5
    chan.start(t_start)
    for t in (t_start + 1e-8 + np.arange(0, 3, 0.05)):
        chan.forward(t)
        chunks = chan.get_new_chunks()[0]
        if len(chunks) > 0:
            print(f"{t}: {chunks}")#, end=", ")

    str(chan.cur_elem)
    str(chan)
    
def perform_random_channel_ops(chan: Channel, random_state: np.random.Generator, max_iterations: int = 1000):
    # issues random operations on the channel
    # random_state: the state to use to sample operations to perform
    # returns: number of actions
    
    nb_actions = {"forward": 0, "stop_receiving": 0, "stop_receiving_missed": 0, "user_unblock": 0, "mux_unblock": 0, "sim_stopped_unblock": 0, "user_unblock_missed": 0, "get_new_chunks": 0, "run_mux_scan": 0, "done": 0}
    i = 0
    while i < max_iterations:
        i += 1
        weights = np.array([0.9, 0.1, 0.1, 0.08, 0.02])
        u = random_state.choice(["forward", "stop_receiving", "user_unblock", "run_mux_scan", "done"], p=weights/weights.sum())
        nb_actions[u] += 1
        if u == "forward":
            chan.forward(random_state.uniform(0.3, 0.6), delta=True)
        elif u == "stop_receiving":
            response = chan.stop_receiving()
            if response == StoppedReceivingResponse.MISSED:
                nb_actions["stop_receiving"] -= 1
                nb_actions["stop_receiving_missed"] += 1
            elif response == StoppedReceivingResponse.ALREADY_STOPPED_RECEIVING:
                nb_actions["stop_receiving"] -= 1
        elif u == "user_unblock":
            missed = not chan.unblock(unblock_duration=random_state.uniform(0.3, 0.6))
            if missed:
                nb_actions["user_unblock_missed"] += 1
                nb_actions["user_unblock"] -= 1
            # else:
            #     print(f"unblocked {read_id}")
        elif u == "run_mux_scan":
            nb_actions["mux_unblock"] += int(isinstance(chan.cur_elem, ChunkedRead))
            
            t_duration = random_state.uniform(0.3, 0.6)
            chan.run_mux_scan(t_duration=t_duration)
            chan.forward(t_duration + 0.01, delta=True)
        
        if random_state.uniform() < 0.2:
            chan.get_new_chunks()
    
    # finish simulation
    nb_actions["sim_stopped_unblock"] += int(isinstance(chan.cur_elem, ChunkedRead))
        
    return nb_actions
        
def test_deterministic(sim_params):
    # test that the channel produces the same results when run twice
    
    def run_channel(seed, chan=None):
        read_pool = ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3)))
        reads_writer = ArrayReadsWriter()
        
        if chan is None:
            chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
            chan.save_elems = True
        else:
            assert not chan.is_running
            assert chan.save_elems
            chan.read_pool = read_pool
            chan.reads_writer = reads_writer
            
        chan.start(0)
        
        nb_actions = perform_random_channel_ops(chan, random_state=np.random.default_rng(seed))
        return (chan.finished_elems + [chan.cur_elem], nb_actions["forward"], chan)
    
    for seed in range(4):
        assert run_channel(seed=seed)[:2] == run_channel(seed=seed)[:2]
    
    # check the same happens when reusing the same channel
    chan_elems, nb_forwards, chan = run_channel(seed=5)
    chan.stop()
    chan_elems2, nb_forwards2 = run_channel(seed=5, chan=chan)[:2]
    chan.stop()
    assert (chan_elems, nb_forwards) == (chan_elems2, nb_forwards2)
    
def test_random_operations(sim_params, channel_write_zero_length_reads):
    # try out a wide range of parameters with start, stop etc and check no assertion errors occur, and plotting
    
    # make all elements roughly the same length, reads have length 8-17 -> take about 1.2s
    sim_params = SimParams(
        gap_samplers={"channel1": ConstantGapsUntilBlocked(short_gap_length=1.2, long_gap_length=1.2, prob_long_gap=0.15, time_until_blocked=np.inf, read_delay=0)},
        bp_per_second=10, chunk_size=4, default_unblock_duration=1.2, seed=0,
    )
    
    plotted_once = False
    for i in range(4):
        all_read_ids = []
        def get_read_and_save_id():
            for (read_id, read_seq) in random_reads_gen(random_state=np.random.default_rng(i)):
                all_read_ids.append(read_id)
                yield read_id, read_seq
                
        read_pool = ReadPoolFromIterable(get_read_and_save_id())
        reads_writer = ArrayReadsWriter()
        chan = Channel("channel1", read_pool, reads_writer, sim_params=sim_params)
        chan.save_elems = True
        
        cum_nb_actions = None
        for j in range(3):
            print(f"i={i}, j={j}", end=", ")
            
            t_start = 2
            chan.start(t_start)
            nb_actions = perform_random_channel_ops(chan, random_state=np.random.default_rng(i + j + 1))
            if cum_nb_actions is None:
                cum_nb_actions = nb_actions
            else:
                cum_nb_actions = {k: cum_nb_actions[k] + nb_actions[k] for k in cum_nb_actions}
            
            if nb_actions["forward"] > 0:
                if not plotted_once:
                    plotted_once = True
                    ax = chan.plot(); plt.close(ax.figure)
                    # plt.show(block=True)
                
            chan.stop()
            
            assert chan.stats.time_active + chan.stats.no_reads_left.time_spent == approx(chan.t - t_start, rel=1e-3)
            
        # check that all reads that were requested by channel were used in this order
        assert [read_id for (read_id, seq, _) in reads_writer.reads] == all_read_ids
        assert sum(["user_unblocked" in x[2] for x in reads_writer.reads]) == cum_nb_actions["user_unblock"]
        assert sum(["stopped_receiving" in x[2] for x in reads_writer.reads]) == cum_nb_actions["stop_receiving"]
        
        #[x[0] for x in reads_writer.reads if "user_unblocked" in x[2]]
        #[x for x in reads_writer.reads if x[0] == "read14"][0]
        