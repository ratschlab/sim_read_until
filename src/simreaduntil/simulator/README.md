# Simulator Structure

It emulates an ONT device controlled by the ReadUntil API directly or via gRPC. It takes full-length reads and terminates them early, as necessary.

## Class Hierarchy

The class hierarchy is as follows
```{bash}
# A -> B means that B inherits from A.

ChannelElement -> ShortGap, LongGap, MuxScan, UnblockDelay, ChunkedRead, NoReadLeftGap

ChunkedRead: an active read that is split into chunks

Channel: repeatedly selects a new channel element

ReadPool -> ReadPoolFromIterable, ReadPoolFromFile, ReadPoolFromIterable

ReadsWriter -> SingleFileReadsWriter, ArrayReadsWriter, RotatingFileReadsWriter

ONTDeviceServicer: gRPC servicer that delegates calls to an underlying ONTSimulator device

ReadUntilDevice: abstract class representing device that allows to start/stop an experiment, block/stop receiving reads, get read chunks
ReadUntilDevice -> ONTSimulator: simulating a given number of channels
               \-> SimulatorGRPCClient: takes server port as argument, connects to a gRPC server that has service ONTDeviceServicer
```

## Developer Notes

We use gRPC because it makes serialization easy, is efficient, allows streaming, and works across languages.
After changing the `.proto` file, you can recompile it with:

```{bash}
# activate python environment with gRPC
# requires gnu sed (brew install gnu-sed) or use: sed -i ''

cd src/simreaduntil/simulator
rm -rf protos_generated
mkdir protos_generated
# replace imports so they work with the project structure ("ont_simulator" must be on the python path)
python -m grpc_tools.protoc -Iprotos/ --python_out=protos_generated/ --pyi_out=protos_generated/ --grpc_python_out=protos_generated/ protos/ont_device.proto && \
    sed -i -E "s%import (.*)_pb2 as%import simreaduntil.simulator.protos_generated.\1_pb2 as%g" protos_generated/ont_device_pb2_grpc.py

# todo: check
cd src && python -m grpc_tools.protoc -Isimreaduntil/simulator/protos/ --python_out=simreaduntil/simulator/protos_generated/ --pyi_out=simreaduntil/simulator/protos_generated/ --grpc_python_out=simreaduntil/simulator/protos_generated/ simreaduntil/simulator/protos/ont_device.proto
```


### Debugging the gRPC server

You can run the tests with `pytest` to see whether it is working. Try to execute the files without pytest to check that the imports work.

You can check that gRPC networking works by looking at the `example_grpc` directory.

Moreover, gRPC logging can be enabled with

```{bash}
export GRPC_VERBOSITY=DEBUG
export GRPC_TRACE=all
```

It may happen that the previous gRPC server instantiation did not terminate properly. In this case, you can kill the process with:

```{bash}
# find relevant python processes via
lsof -i tcp | grep Python

# if you know the port, e.g. 10871
netstat -an | grep 10871
sudo lsof -i :10871
#sudo kill -9 <PID>
```