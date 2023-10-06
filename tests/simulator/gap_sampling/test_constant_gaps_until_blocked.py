import numpy as np

from simreaduntil.simulator.channel_stats import ChannelStats
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.gap_sampling.gap_sampling import GapSampler


def test_constant_gap_sampler():
    gap_sampler = ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0)
    random_state = np.random.default_rng(1)
    channel_stats = ChannelStats()
    for i in range(10):
        gap_sampler.sample_next_gap(channel_stats, random_state=random_state)
        channel_stats.reads.add_full(time=2, nb_bps=10, stopped_receiving=False)
    assert gap_sampler.sample_read_start_delay(channel_stats, random_state=random_state) == 0

    gap_sampler = ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=1.0, time_until_blocked=24, read_delay=0)
    channel_stats = ChannelStats()
    gap_type, _ = gap_sampler.sample_next_gap(channel_stats, random_state=random_state)
    assert gap_type == GapSampler.GapType.Long
    channel_stats.reads.add_full(time=25, nb_bps=5, stopped_receiving=False)
    gap_type, _ = gap_sampler.sample_next_gap(channel_stats, random_state=random_state)
    assert gap_type == GapSampler.GapType.Broken
    