from pathlib import Path
import grpc
import numpy as np
import pytest

from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger
from simreaduntil.simulator.gap_sampling.constant_gaps_until_blocked import ConstantGapsUntilBlocked
from simreaduntil.simulator.readpool import ReadPoolFromIterable
from simreaduntil.simulator.readswriter import ArrayReadsWriter
from simreaduntil.simulator.simulator import ONTSimulator
from simreaduntil.simulator.simulator_params import SimParams
import simreaduntil.simulator.simulator_server as sim_server
from simreaduntil.simulator.simulator_server import launchable_device_grpc_server, manage_grpc_server, port_is_available, print_gen_exceptions, print_nongen_exceptions
from simreaduntil.simulator.protos_generated import ont_device_pb2_grpc, ont_device_pb2
from simreaduntil.usecase_helpers.utils import random_reads_gen


# capsys to capture stdout
def test_print_nongen_exceptions(capsys):
    
    def f_normal():
        a = 2
        return a * 4
    
    def f_error():
        a = 2
        raise ValueError("error:ferror")
        return a
    
    def gen_normal():
        yield from range(3)
        
    def gen1():
        yield from range(3)
        raise ValueError("error:gen1")
    
    # final yield
    def gen2():
        yield from range(3)
        raise ValueError("error:gen2")
        yield 3
        
        
    # calling it with generator instead of normal function and vice versa
    with pytest.raises(AssertionError, match="must be a generator"):
        print_gen_exceptions(f_error)
    with pytest.raises(AssertionError, match="must not be a generator"):
        print_nongen_exceptions(gen1)
        
    # exception will only be printed, not reraised
    h = print_nongen_exceptions(f_error)
    h()
    
    h_gen = print_gen_exceptions(gen1)
    list(h_gen())
    
    # also raise exceptions that occur
    sim_server.RAISE_GRPC_SERVER_EXCEPTIONS = True
    
    h = print_nongen_exceptions(f_normal)
    h()
    
    h = print_nongen_exceptions(f_error)
    with pytest.raises(ValueError, match="error:ferror"):
        h()
        
    h_gen = print_gen_exceptions(gen_normal)
    list(h_gen())
    
    h_gen = print_gen_exceptions(gen1)
    with pytest.raises(ValueError, match="error:gen1"):
        list(h_gen())
        
    h_gen = print_gen_exceptions(gen2)
    with pytest.raises(ValueError, match="error:gen2"):
        list(h_gen())
    
        
    capsys.readouterr() # otherwise stderr is printed

def test_launchable_device_grpc_server():
    sim_params = SimParams(
        gap_samplers={f"channel_{i}": ConstantGapsUntilBlocked(short_gap_length=0.4, long_gap_length=10.1, prob_long_gap=0, time_until_blocked=np.inf, read_delay=0) for i in range(2)},
    )
    simulator = ONTSimulator(
        read_pool=ReadPoolFromIterable(random_reads_gen(random_state=np.random.default_rng(3), length_range=(10, 50))), 
        reads_writer=ArrayReadsWriter(),
        sim_params=sim_params
    )
    
    # import os; os.environ["GRPC_VERBOSITY"] = "DEBUG"; os.environ["GRPC_TRACE"] = "http"
    
    port, server, unique_id = launchable_device_grpc_server(simulator)
    assert port != 0
    assert not port_is_available(port)
    
    with manage_grpc_server(server):
        with grpc.insecure_channel(f'localhost:{port}') as channel:
            stub = ont_device_pb2_grpc.ONTDeviceStub(channel)
            
            assert stub.GetServerInfo(ont_device_pb2.EmptyRequest()).unique_id == unique_id, "mismatching unique_ids, probably connected to wrong server"
            assert not stub.IsRunning(ont_device_pb2.EmptyRequest()).value
            
            stub.GetActionResults(ont_device_pb2.ActionResultsRequest(clear=True))
            
            assert stub.StartSim(ont_device_pb2.StartRequest(acceleration_factor=2)).value
            
            # unblocking inexistent read
            assert not stub.PerformActions(ont_device_pb2.ReadActionsRequest(actions=[
                ont_device_pb2.ReadActionsRequest.Action(channel=2, read_id="inexistent", unblock=ont_device_pb2.ReadActionsRequest.Action.UnblockAction(unblock_duration=0.2))
            ])).succeeded[0]
            
            assert not stub.PerformActions(ont_device_pb2.ReadActionsRequest(actions=[
                ont_device_pb2.ReadActionsRequest.Action(channel=1, read_id="inexistent", stop_further_data=ont_device_pb2.ReadActionsRequest.Action.StopReceivingAction()),
            ])).succeeded[0]
            
            assert stub.StopSim(ont_device_pb2.EmptyRequest()).value
            
    # assert port_is_available(port) # port not yet released
    
    
# add_comprehensive_stream_handler_to_logger()
# test_launch_device_grpc_server()

