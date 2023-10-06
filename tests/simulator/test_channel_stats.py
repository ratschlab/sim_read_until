
import pytest
from pytest import approx
from simreaduntil.simulator.channel_stats import ChannelStats, ElementStats, ReadElementStats, combine_stats, plot_channel_stats

dict_different_keys = pytest.helpers.dict_different_keys


def test_ChannelStats():
    channel_stats = ChannelStats(n_channels=1)
    channel_stats.check_consistent()
    assert channel_stats.n_channels_running == 0
    
    channel_stats.short_gaps.cur_number = 1
    channel_stats.check_consistent()
    assert channel_stats.n_channels_running == 1

    # test printing
    print(channel_stats)
    channel_stats.short_gaps.time_spent = 2.1
    print(channel_stats)
    
def test_channelstat_methods():
    channel_stats = ChannelStats(
        n_channels=14,
        short_gaps=ElementStats(cur_number=2, finished_number=11, time_spent=4.2), # 11 >= 8+2
        long_gaps=ElementStats(cur_number=1, finished_number=1, time_spent=2.3),
        unblock_delays=ElementStats(cur_number=2, finished_number=1, time_spent=1.2),
        mux_scans=ElementStats(cur_number=1, finished_number=2),
        reads=ReadElementStats(
            cur_number=2, finished_number=8, time_spent=5.6,
            cur_number_stop_receiving=1,
            fin_number_rejected=3, number_rejected_missed=1, 
            fin_number_stop_receiving=1, number_stop_receiving_missed=2,
            number_bps_read=240, number_bps_rejected=120, number_bps_requested=130,
        ),
        no_reads_left=ElementStats(cur_number=2, time_spent=11.1),
    )
    
    channel_stats.check_consistent()
    assert channel_stats.n_channels_running == 10
    time_active = 4.2 + 2.3 + 1.2 + 5.6
    assert channel_stats.time_active == time_active
    
    assert channel_stats.number_full_reads() == 8 - 3
    assert channel_stats.fraction_reading_time() == approx(5.6 / time_active)
    assert channel_stats.fraction_short_gaps() == approx(4.2 / time_active)
    assert channel_stats.fraction_unblocking() == approx(1.2 / time_active)
    assert channel_stats.fraction_long_gaps() == approx(2.3 / time_active)
    
    # check_consistent
    channel_stats = ChannelStats(n_channels=1)
    with pytest.raises(AssertionError):
        channel_stats.short_gaps.cur_number = 1
        channel_stats.reads.number_bps_read = 120
        channel_stats.reads.number_bps_requested = 130 # more processed than read
        channel_stats.check_consistent()
    
def test_start_add_finish_element():
    channel_stats = ChannelStats(n_channels=1)
    
    channel_stats.mux_scans.start(); channel_stats.check_consistent()
    channel_stats.mux_scans.add_time(0.5); channel_stats.check_consistent()
    channel_stats.mux_scans.finish(); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.short_gaps.start(); channel_stats.check_consistent()
    channel_stats.short_gaps.add_time(2.1); channel_stats.check_consistent()
    channel_stats.short_gaps.add_time(2.1); channel_stats.check_consistent()
    channel_stats.short_gaps.finish(); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.reads.start(); channel_stats.check_consistent()
    channel_stats.reads.add_time(1.4, nb_new_bps=10); channel_stats.check_consistent()
    channel_stats.reads.finish(stopped_receiving=False, nb_bps_rejected=3); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.unblock_delays.start(); channel_stats.check_consistent()
    channel_stats.unblock_delays.add_time(0.3); channel_stats.check_consistent()
    channel_stats.unblock_delays.finish(); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.long_gaps.start(); channel_stats.check_consistent()
    channel_stats.long_gaps.add_time(1.0); channel_stats.check_consistent()
    channel_stats.long_gaps.add_time(1.2); channel_stats.check_consistent()
    channel_stats.long_gaps.finish(); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.reads.start(); channel_stats.check_consistent()
    channel_stats.reads.add_time(1.9, nb_new_bps=12); channel_stats.check_consistent()
    channel_stats.reads.cur_number_stop_receiving += 1
    channel_stats.reads.finish(stopped_receiving=True, nb_bps_rejected=0); channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    channel_stats.no_reads_left.start(); channel_stats.check_consistent()
    channel_stats.no_reads_left.add_time(0.5); channel_stats.check_consistent()
    channel_stats.no_reads_left.finish()
    channel_stats.check_consistent(n_channels_running=0)
    print(str(channel_stats))
    
    expected_channel_stats = ChannelStats(
        n_channels=1,
        short_gaps=ElementStats(finished_number=1, time_spent=2.1 + 2.1),
        long_gaps=ElementStats(finished_number=1, time_spent=1.0 + 1.2),
        unblock_delays=ElementStats(finished_number=1, time_spent=0.3),
        mux_scans=ElementStats(finished_number=1, time_spent=0.5),
        no_reads_left=ElementStats(finished_number=1, time_spent=0.5),
        reads=ReadElementStats(finished_number=2, time_spent=1.4 + 1.9, fin_number_rejected=1, fin_number_stop_receiving=1, number_bps_read=10 + 12, number_bps_rejected=3),
    )
    assert channel_stats.n_channels_running == 0
    assert channel_stats.time_active == 0.5 + 2.1 + 2.1 + 1.4 + 0.3 + 1.0 + 1.2 + 1.9
    assert channel_stats == expected_channel_stats
    
    # check channel_broken
    channel_stats = ChannelStats(
        n_channels=5,
        short_gaps=ElementStats(cur_number=3, finished_number=2, time_spent=2.1 + 2.1),
        channel_broken=ElementStats(cur_number=1, finished_number=2, time_spent=1.0 + 1.2),
    )
    channel_stats.check_consistent()
    
def test_combine_stats():
    channel_stats1 = ChannelStats(
        n_channels=3,
        short_gaps=ElementStats(cur_number=1, finished_number=2, time_spent=2.1),
        reads=ReadElementStats(cur_number=2),
    )
    channel_stats1.check_consistent()
    
    channel_stats2 = ChannelStats(
        n_channels=6,
        short_gaps=ElementStats(cur_number=2, finished_number=3, time_spent=4.4),
        reads=ReadElementStats(cur_number=3),
    )
    channel_stats2.check_consistent()
    
    combined_stats = combine_stats((channel_stats1, channel_stats2))
    assert combined_stats.n_channels_running == 3 + 5
    assert combined_stats.n_channels == 3 + 6
    assert combined_stats.short_gaps.cur_number == 3
    assert combined_stats.short_gaps.time_spent == 6.5
    
def test_plot_channel_stats():
    channel_stats = [
        ChannelStats(
            n_channels=1,
            short_gaps=ElementStats(finished_number=450, time_spent=30),
            long_gaps=ElementStats(finished_number=20, time_spent=45),
            unblock_delays=ElementStats(finished_number=30, time_spent=20),
            mux_scans=ElementStats(finished_number=8, time_spent=25),
            channel_broken=ElementStats(finished_number=5, time_spent=20),
            reads=ReadElementStats(finished_number=420, time_spent=100, number_bps_requested=1200, number_bps_read=1800, number_bps_rejected=3200, fin_number_stop_receiving=12, fin_number_rejected=33, number_rejected_missed=11, number_stop_receiving_missed=9),
            no_reads_left=ElementStats(finished_number=3, time_spent=15),
        ),
        ChannelStats(
            n_channels=1,
            short_gaps=ElementStats(finished_number=250, time_spent=20),
            long_gaps=ElementStats(finished_number=40, time_spent=25),
            unblock_delays=ElementStats(finished_number=10, time_spent=40),
            mux_scans=ElementStats(finished_number=18, time_spent=30),
            channel_broken=ElementStats(finished_number=10, time_spent=10),
            reads=ReadElementStats(finished_number=240, time_spent=100, number_bps_requested=700, number_bps_read=1000, number_bps_rejected=500, fin_number_stop_receiving=5, fin_number_rejected=25, number_rejected_missed=15, number_stop_receiving_missed=20),
            no_reads_left=ElementStats(finished_number=7, time_spent=35),
        )
    ]
    plot_channel_stats(channel_stats)
