"""
gRPC Client to connect to a gRPC server exposing an ONTSimulator
"""

import grpc
from typing import Any, List, Optional, Tuple

from simreaduntil.simulator.channel import StoppedReceivingResponse, UnblockResponse
from simreaduntil.simulator.simulator import ActionType
from simreaduntil.simulator.simulator import ReadUntilDevice
from simreaduntil.simulator.protos_generated import ont_device_pb2_grpc, ont_device_pb2
from simreaduntil.shared_utils.utils import setup_logger_simple

logger = setup_logger_simple(__name__)
"""module logger"""
 
class DeviceGRPCClient(ReadUntilDevice):
    """
    Create a virtual device by connecting to a gRPC server run by a SimulatorGRPCServer
    
    Before calling the device functions, you first have to connect. This class supports the with statement.
    For details about these functions, look at the underlying device's documentation (that you connect to).
    This wraps the gRPC stub class.
    """
    def __init__(self, port, host="localhost"):
        assert port > 0
        self._port : int = port
        self._host : str = host
        
        self._stub : Optional[ont_device_pb2_grpc.ONTDeviceStub] = None
        self._channel : Optional[grpc.Channel] = None
    
    ##### Connection-specific methods
    
    def connect(self):
        logger.info("Connecting to simulator server")
        self._channel = grpc.insecure_channel(f'{self._host}:{self._port}')
        self._stub = ont_device_pb2_grpc.ONTDeviceStub(self._channel)
        logger.info("Connected to simulator server")
        
    @property
    def is_connected(self):
        return self._stub is not None
    
    def _check_connected(self):
        assert self.is_connected, "Simulator client not connected to server"
    
    def disconnect(self):
        logger.info("Disconnecting from simulator server")
        self._stub = None
        self._channel.close()
        self._channel = None
        logger.info("Disconnected from simulator server")
        
    # with statement
    def __enter__(self):
        self.connect()
        self.wait_for_ready()
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.disconnect()
    
    def wait_for_ready(self):
        """Wait until the gRPC server is ready"""
        self._check_connected()
        # statement only finishes once the server is running, response value does not matter here
        self._stub.IsRunning(ont_device_pb2.EmptyRequest(), wait_for_ready=True)
    
    ##### Simulator-specific methods
     
    @property
    def n_channels(self):
        """
        Number of channels
        """
        self._check_connected()
        return self._stub.GetDeviceInfo(ont_device_pb2.EmptyRequest()).n_channels
    
    def device_info(self) -> str:
        """
        print information about the device
        """
        self._check_connected()
        return self._stub.GetDeviceInfo(ont_device_pb2.EmptyRequest()).info
    
    @property
    def unique_id(self) -> str:
        """
        Unique ID of the server, to check that one is connected to the correct server
        """
        return self._stub.GetServerInfo(ont_device_pb2.EmptyRequest()).unique_id
        
    def start(self, acceleration_factor: float = 1.0):
        """
        Start the sequencing.
        """
        self._check_connected()
        return self._stub.StartSim(ont_device_pb2.StartRequest(acceleration_factor=acceleration_factor)).value
    
    def stop(self):
        """
        Stop the sequencing. Write all partial reads to a file (flush).
        """
        self._check_connected()
        return self._stub.StopSim(ont_device_pb2.EmptyRequest()).value
    
    def run_mux_scan(self, t_duration: float):
        """
        Run mux scan
        """
        self._check_connected()
        return self._stub.RunMuxScan(ont_device_pb2.RunMuxScanRequest()).nb_reads_rejected
    
    @property
    def is_running(self):
        """
        Whether the device is sequencing
        """
        self._check_connected()
        return self._stub.IsRunning(ont_device_pb2.EmptyRequest()).value
    
    def get_basecalled_read_chunks(self, batch_size=None, channel_subset=None):
        """
        Get available read chunks from the selected channels, from at most 'batch_size' channels
        """
        self._check_connected()
        if channel_subset is None:
            channels = None
        else:
            channels = ont_device_pb2.BasecalledChunksRequest.Channels(value=channel_subset)
        for chunk in self._stub.GetBasecalledChunks(ont_device_pb2.BasecalledChunksRequest(batch_size=batch_size, channels=channels)):
            yield (chunk.channel, chunk.read_id, chunk.seq, chunk.quality_seq, chunk.estimated_ref_len_so_far)
            
    def get_action_results(self, clear=True) -> List[Tuple[Any, float, int, str, Any]]:
        """
        Get action results
        """
        for action_response in self._stub.GetActionResults(ont_device_pb2.ActionResultsRequest(clear=clear)).actions:
            action_type = ActionType(action_response.action_type)
            action_result = (StoppedReceivingResponse if action_type == ActionType.StopReceiving else UnblockResponse)(action_response.result)
            yield (action_response.read_id, action_response.time, action_response.channel, action_type, action_result)
    
    def unblock_read(self, read_channel, read_id, unblock_duration=None):
        """
        Unblock read_id on channel; returns whether the action was performed (not performed if the read was already over)
        """
        self._check_connected()
        return self._stub.PerformActions(ont_device_pb2.ReadActionsRequest(actions=[
            ont_device_pb2.ReadActionsRequest.Action(channel=read_channel, read_id=read_id, unblock=ont_device_pb2.ReadActionsRequest.Action.UnblockAction(unblock_duration=unblock_duration))    
        ])).succeeded[0]
        
    def stop_receiving_read(self, read_channel, read_id):
        """
        Stop receiving read_id on channel; returns whether the action was performed (not performed if the read was already over)
        """
        self._check_connected()
        return self._stub.PerformActions(ont_device_pb2.ReadActionsRequest(actions=[
            ont_device_pb2.ReadActionsRequest.Action(channel=read_channel, read_id=read_id, stop_further_data=ont_device_pb2.ReadActionsRequest.Action.StopReceivingAction()),
        ])).succeeded[0]
    
    @property
    def mk_run_dir(self):
        """Directory where reads are written to"""
        self._check_connected()
        return self._stub.GetMKRunDir(ont_device_pb2.EmptyRequest()).mk_run_dir
    