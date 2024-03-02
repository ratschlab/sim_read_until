
import argparse
import logging
import signal
import time
import grpc

import numpy as np
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, print_logging_levels, setup_logger_simple

from simreaduntil.shared_utils.utils import print_args
from simreaduntil.simulator.simulator_client import DeviceGRPCClient
from simreaduntil.shared_utils.utils import set_signal_handler

logger = setup_logger_simple(__name__)
"""module logger"""

def parse_args():
    parser = argparse.ArgumentParser(description="Simple gRPC ReadUntil client that randomly decides to stop receiving, reject reads or take no action")
    parser.add_argument("port", type=int, help="Port to connect to")
    
    args = parser.parse_args()
    print_args(args, logger=logger)
    return args

def main():
    log_level = logging.INFO
    logging.getLogger(None).setLevel(log_level)
    add_comprehensive_stream_handler_to_logger(None, level=log_level)
    # print_logging_levels()
    
    args = parse_args()
    
    port = args.port
    assert port > 0, f"port {port} must be > 0"
    
    rng = np.random.default_rng(42)
    
    with DeviceGRPCClient(port) as client:
        assert client.is_connected, "client not connected"
        
        if not client.start(acceleration_factor=1):
            logger.warning("Simulator was already running")
        assert client.is_running, "simulator not running"
        
        logger.info(f"Simulator writes reads to directory: {client.mk_run_dir}")
        logger.info(f"Simulator has the following properties: {client.device_info()}")
        
        num_batches = 0
        num_chunks = 0
        
        def stop_client(*args, **kwargs):
            try:
                if client.stop():
                    logger.info("Stopped simulation")
            except grpc.RpcError as e:
                pass
        
        with set_signal_handler(signal_type=signal.SIGINT, handler=stop_client): # catch keyboard interrupt (Ctrl+C)
            try:
                with logging_redirect_tqdm():
                    while client.is_running:
                        num_batches += 1
                        for (channel, read_id, seq, quality_seq, estimated_ref_len_so_far) in tqdm(client.get_basecalled_read_chunks(), desc=f"Processing chunks in batch {num_batches}"):
                            num_chunks += 1
                            logger.debug(f"Read chunk: channel={channel}, read_id={read_id}, seq={seq[:20]}..., quality_seq={quality_seq}, estimated_ref_len_so_far={estimated_ref_len_so_far}")
                            u = rng.uniform()
                            if u < 0.2:
                                logger.debug(f"Rejecting read '{read_id}'")
                                client.unblock_read(channel, read_id)
                            elif u < 0.4:
                                logger.debug(f"Stop receiving read '{read_id}'")
                                client.stop_receiving_read(channel, read_id)
                            else:
                                # no action
                                pass
                            # time.sleep(0.05)
                        time.sleep(0.2) # throttle
            except grpc.RpcError as e:
                logger.error(f"Caught gRPC error: {e}")
            
            
    logger.info(f"Done. Received {num_chunks} chunks from {num_batches} batches")
        
        
if __name__ == "__main__":
    main()