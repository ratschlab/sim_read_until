import numpy as np
import pytest
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.readpool import ReadPoolFromIterable
from simreaduntil.simulator.readswriter import ArrayReadsWriter
from simreaduntil.simulator.simulator import ONTSimulator
from simreaduntil.simulator.simulator_client import DeviceGRPCClient
from simreaduntil.simulator.simulator_params import SimParams
from simreaduntil.usecase_helpers.utils import random_reads_gen


from simreaduntil.simulator.simulator_server import launchable_device_grpc_server, manage_grpc_server

check_simulator_actions_agree_with_reads_writer = pytest.helpers.check_simulator_actions_agree_with_reads_writer
perform_random_sim_ops = pytest.helpers.perform_random_sim_ops

            
def test_grpc_client(channel_write_zero_length_reads):
    reads_writer = ArrayReadsWriter()
    reads_writer.output_dir = "dummy_dir" # patch attribute
    
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
    )
    simulator = ONTSimulator(
        read_pool=ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3), length_range=(10, 50))), 
        reads_writer=reads_writer,
        sim_params=sim_params,
    )
    
    port, server, unique_id = launchable_device_grpc_server(simulator)
    assert port != 0
        
    acceleration_factor = 5
    seed = 2
    
    with manage_grpc_server(server):
        with DeviceGRPCClient(port) as client:
            assert client.is_connected
            
            assert client.unique_id == unique_id, f"mismatching unique_ids, probably connected to an existing server: {client.unique_id} != {unique_id}"
            
            assert str(client.mk_run_dir) == "dummy_dir"
            assert not client.is_running
            
            assert client.start(acceleration_factor=acceleration_factor)
            assert not client.start(), "already started"
            
            assert client.n_channels == simulator.n_channels
            assert isinstance(client.device_info(), str)
            
            nb_actions = perform_random_sim_ops(client, random_state=np.random.default_rng(seed), acceleration_factor=acceleration_factor, async_mode=True, nb_iterations=100)
            
            assert client.stop()
            assert not client.stop(), "already stopped"
            
            check_simulator_actions_agree_with_reads_writer(simulator, nb_actions)
            
            assert client.is_connected
            
        assert not client.is_connected
        
        client = DeviceGRPCClient(port)
        assert not client.is_connected
