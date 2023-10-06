# Example gRPC

This can be used to check that the network settings work. Code is taken from grpc routeguide example.

Run as follows in folder containing proto file

```{bash}
rm -rf example_grpc/generated
mkdir example_grpc/generated
(cd example_grpc/generated && python -m grpc_tools.protoc -I../ --python_out=. --pyi_out=. --grpc_python_out=. ../routeguide.proto)
(cd example_grpc/generated && sed -i '' -E "s%import (.*)_pb2 as%import example_grpc.generated.\1_pb2 as%g" routeguide_pb2_grpc.py)
# activate python environment with grpc
python -m example_grpc.example_server &; python -m example_grpc.example_client
```

<!-- The python path is incorrect for the following, so example_grpc is not found
python example_grpc/example_server.py &; python example_grpc/example_client.py -->