"""
Generate an example plot for the simulator with 2 channels (included in the paper)
"""

import numpy as np
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.readpool import ReadPoolFromIterable
from simreaduntil.simulator.readswriter import ArrayReadsWriter
from simreaduntil.simulator.simulator import ONTSimulator
from simreaduntil.simulator.simulator_params import SimParams

sim_params = SimParams(
    gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=1.5, prob_long_gap=0.5, time_until_blocked=8.1, read_delay=0) for i in range(2)},
    bp_per_second=10, min_chunk_size=4, default_unblock_duration=0.8, seed=0,
)

rng = np.random.default_rng(0)
def reads_gen():
    for i in range(8):
        l = rng.integers(1, 3)*10
        yield f"read{i}", "A" * l

read_pool = ReadPoolFromIterable(reads_gen())
simulator = ONTSimulator(
    read_pool=read_pool,
    reads_writer=ArrayReadsWriter(),
    sim_params = sim_params,
    output_dir="<unavailable>",
)
simulator.save_elems = True

simulator.sync_start(0)
simulator.sync_forward(2)
simulator._channels[1].unblock()
simulator.sync_forward(5.5)
# simulator._channels[0].cur_elem.
simulator._channels[0].unblock()
simulator.sync_forward(8)
simulator.run_mux_scan(2, is_sync=True)
simulator.sync_forward(13)

ax = simulator.plot_channels()#; import matplotlib.pyplot as plt; plt.show()
ax.figure.tight_layout()
ax.set_ylim([-0.2, 1.2])
ax.autoscale()

simulator.sync_stop()
ax.figure.savefig("simulator_example.png", dpi=300)