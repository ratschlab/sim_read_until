"""
gRPC server exposing an underlying ReadUntilDevice

API calls: start/stop, unblock/stop_receiving, mk_run_dir, is_running, basecalled_chunks
"""

from contextlib import contextmanager
from pathlib import Path
import argparse
import shutil
import tempfile
import time
import logging
from typing import Tuple
import unittest
from unittest import mock
import grpc
from concurrent import futures
from simreaduntil.simulator.channel import StoppedReceivingResponse

from simreaduntil.simulator.simulator import ReadUntilDevice, ONTSimulator
from simreaduntil.simulator.protos_generated import ont_device_pb2_grpc, ont_device_pb2
from simreaduntil.shared_utils.utils import setup_logger_simple


logger = setup_logger_simple(__name__)
"""module logger"""

RAISE_GRPC_SERVER_EXCEPTIONS = False
"""Whether to only print gRPC server exceptions or also raise them"""

def print_nongen_exceptions(f):
    """
    gRPC catches exceptions, so this function makes them visible before they are caught and then raises them again if RAISE_GRPC_SERVER_EXCEPTIONS is set
    
    Does not work when f is a generator function because it may not have finished when it returns (paused at the yield statement).
    See print_generator_exceptions
    """
    import inspect
    assert not inspect.isgeneratorfunction(f), "must not be a generator"
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            import traceback
            traceback.print_exc()
            if RAISE_GRPC_SERVER_EXCEPTIONS:
                raise
    inner.__name__ = f.__name__
    return inner

def print_gen_exceptions(f):
    """
    See print_nongen_exceptions, works with generators
    """
    import inspect
    assert inspect.isgeneratorfunction(f), "must be a generator"
    def inner(*args, **kwargs):
        try:
            yield from f(*args, **kwargs)
        except Exception as e:
            import traceback
            traceback.print_exc()
            yield # needed for exception to be handled, see https://mail.python.org/pipermail/python-list/2003-June/216000.html
            if RAISE_GRPC_SERVER_EXCEPTIONS:
                raise
    inner.__name__ = f.__name__
    return inner

# See example at https://github.com/grpc/grpc/blob/v1.50.0/examples/protos/route_guide.proto and https://github.com/grpc/grpc/blob/v1.50.0/examples/python/route_guide/route_guide_server.py
class ONTDeviceServicer(ont_device_pb2_grpc.ONTDeviceServicer):
    """
    gRPC server that exposes the simulator via a gRPC server
    
    Use launch_device_grpc_server() to launch such a server
    
    Args:
        device: ReadUntilDevice to expose via gRPC
        unique_id: unique id, this can be set to the start time to check that one is interacting with the correct server
    """
    def __init__(self, device : ReadUntilDevice, unique_id: str="notset"):
        super().__init__()
        self.device = device
        self.unique_id = unique_id
    
    @print_nongen_exceptions
    def GetServerInfo(self, request, context):
        return ont_device_pb2.ServerInfoResponse(unique_id=self.unique_id)
    
    @print_nongen_exceptions
    def GetMKRunDir(self, request, context):
        # need str(), otherwise _InactiveRpcError
        return ont_device_pb2.MKRunDirResponse(mk_run_dir=str(self.device.mk_run_dir))
            
    @print_nongen_exceptions
    def PerformActions(self, request, context):
        res = []
        for action_desc in request.actions:
            channel, read_id = action_desc.channel, action_desc.read_id
            if action_desc.WhichOneof("action") == "unblock":
                unblock_duration = None if action_desc.unblock.unblock_duration < 0 else action_desc.unblock.unblock_duration
                res.append(self.device.unblock_read(channel, read_id=read_id, unblock_duration=unblock_duration))
            else:
                res.append(self.device.stop_receiving_read(channel, read_id=read_id)) #todo2: current conversion from enum 0,1,2 to bool is not ideal
        return ont_device_pb2.ActionResultImmediateResponse(succeeded=res)
    
    @print_gen_exceptions
    def GetActionResults(self, request, context):
        for (read_id, time, channel, action_type, action_result) in self.device.get_action_results(clear=request.clear):
            yield ont_device_pb2.ActionResultResponse(read_id=read_id, time=time, channel=channel, action_type=int(action_type), action_result=int(action_result))
        
    @print_gen_exceptions
    def GetBasecalledChunks(self, request, context):
        # channel_subset=None on request side means that field was not set
        channel_subset = request.channels.value if request.HasField("channels") else None
        batch_size = request.batch_size if request.batch_size > 0 else None
        for (channel, read_id, seq, quality_seq, estimated_ref_len_so_far) in self.device.get_basecalled_read_chunks(batch_size=batch_size, channel_subset=channel_subset):
            yield ont_device_pb2.BasecalledReadChunkResponse(channel=channel, read_id=read_id, seq=seq, quality_seq=quality_seq, estimated_ref_len_so_far=estimated_ref_len_so_far)

    @print_nongen_exceptions
    def StartSim(self, request, context):
        """
        Start simulation, don't raise an error if the device is not running as this can easily be recovered (no need to crash)
        
        Returns: whether it succeeded (i.e. if simulation was not running)
        """
        acceleration_factor = request.acceleration_factor if request.acceleration_factor <= 0 else 1.0
        return ont_device_pb2.BoolResponse(value=self.device.start(acceleration_factor=acceleration_factor, update_method=request.update_method, log_interval=request.log_interval, stop_if_no_reads=request.stop_if_no_reads))
    
    @print_nongen_exceptions
    # stop simulation, returns whether it succeeded (i.e. if simulation was running)
    def StopSim(self, request, context):
        return ont_device_pb2.BoolResponse(value=self.device.stop())
        
    @print_nongen_exceptions
    def RunMuxScan(self, request, context):
        assert request.HasField("t_duration"), "t_duration must be set"
        return ont_device_pb2.MuxScanStartedInfo(value=self.device.run_mux_scan(t_duration=request.t_duration))
    
    @print_nongen_exceptions
    # whether simulation is running
    def IsRunning(self, request, context):
        return ont_device_pb2.BoolResponse(value=self.device.is_running)
    
    @print_nongen_exceptions
    # get parameters of the simulation
    def GetDeviceInfo(self, request, context):
        if isinstance(self.device, ONTSimulator):
            info = self.device.device_info(sim_params=True, channel_states=False)
        else:
            info = self.device.device_info()
            
        return ont_device_pb2.DeviceInfoResponse(
            # assumes device_info to take these parameters
            info=info,
            n_channels=self.device.n_channels,
        )

def port_is_available(port) -> bool:
    """
    Whether a port is available to be used with the grpc server
    """
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        result = sock.connect_ex(('127.0.0.1', port))
    return result != 0

def launchable_device_grpc_server(device: ReadUntilDevice, port: int=0) -> Tuple[int, grpc.Server, str]:
    """
    Start a grpc server for interacting with the given device on port; adds support for with statement
    
    Args:
        device: underlying device
        port: port where to launch, assigns an available port if set to 0
        
    Returns:
        A tuple containing (port, server, unique_id) where
        - port is the actual port (useful if port set to 0 in argument),
        - server must be started,
        - unique_id that the server will respond with GetUniqueId(); to check that the correct server is being used
    Raises:
        AssertionError if port is already in use
    
    """
    # grpc also starts server even if it is in use! -> may get responses from wrong server, so check for it
    assert port_is_available(port), f"Port {port} is already in use"
    
    # launch server
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=1)) # only use 1 worker because ONTSimulator is not thread-safe (no longer true, need to test his)
    def get_unique_id() -> str:
        import datetime
        import uuid
        return str(f"Started {datetime.datetime.now()}, uuid={uuid.uuid1()}")
    
    unique_id = get_unique_id()
    ont_device_pb2_grpc.add_ONTDeviceServicer_to_server(ONTDeviceServicer(device, unique_id=unique_id), server)
    
    # port returned if port is 0
    port = server.add_insecure_port(f"[::]:{port}")
    logger.info(f"Launched server on port {port}")
    
    assert not port_is_available(port), f"gRPC server should be running on port {port}"
    
    return port, server, unique_id

@contextmanager
def manage_grpc_server(server: grpc.server, grace=0.5):
    """Provides a with context around the gRPC server (not available)"""
    try:
        logger.info("Starting server")
        server.start()
        yield server
    finally:
        logger.info(f"Stopping server within {grace}s")
        server.stop(grace=grace)
