from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Mapping as _Mapping, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class EmptyRequest(_message.Message):
    __slots__ = ()
    def __init__(self) -> None: ...

class EmptyResponse(_message.Message):
    __slots__ = ()
    def __init__(self) -> None: ...

class ServerInfoResponse(_message.Message):
    __slots__ = ("unique_id",)
    UNIQUE_ID_FIELD_NUMBER: _ClassVar[int]
    unique_id: str
    def __init__(self, unique_id: _Optional[str] = ...) -> None: ...

class MKRunDirResponse(_message.Message):
    __slots__ = ("mk_run_dir",)
    MK_RUN_DIR_FIELD_NUMBER: _ClassVar[int]
    mk_run_dir: str
    def __init__(self, mk_run_dir: _Optional[str] = ...) -> None: ...

class ReadActionsRequest(_message.Message):
    __slots__ = ("actions",)
    class Action(_message.Message):
        __slots__ = ("channel", "read_id", "unblock", "stop_further_data")
        class StopReceivingAction(_message.Message):
            __slots__ = ()
            def __init__(self) -> None: ...
        class UnblockAction(_message.Message):
            __slots__ = ("unblock_duration",)
            UNBLOCK_DURATION_FIELD_NUMBER: _ClassVar[int]
            unblock_duration: float
            def __init__(self, unblock_duration: _Optional[float] = ...) -> None: ...
        CHANNEL_FIELD_NUMBER: _ClassVar[int]
        READ_ID_FIELD_NUMBER: _ClassVar[int]
        UNBLOCK_FIELD_NUMBER: _ClassVar[int]
        STOP_FURTHER_DATA_FIELD_NUMBER: _ClassVar[int]
        channel: int
        read_id: str
        unblock: ReadActionsRequest.Action.UnblockAction
        stop_further_data: ReadActionsRequest.Action.StopReceivingAction
        def __init__(self, channel: _Optional[int] = ..., read_id: _Optional[str] = ..., unblock: _Optional[_Union[ReadActionsRequest.Action.UnblockAction, _Mapping]] = ..., stop_further_data: _Optional[_Union[ReadActionsRequest.Action.StopReceivingAction, _Mapping]] = ...) -> None: ...
    ACTIONS_FIELD_NUMBER: _ClassVar[int]
    actions: _containers.RepeatedCompositeFieldContainer[ReadActionsRequest.Action]
    def __init__(self, actions: _Optional[_Iterable[_Union[ReadActionsRequest.Action, _Mapping]]] = ...) -> None: ...

class ActionResultsRequest(_message.Message):
    __slots__ = ("clear",)
    CLEAR_FIELD_NUMBER: _ClassVar[int]
    clear: bool
    def __init__(self, clear: bool = ...) -> None: ...

class ActionResultResponse(_message.Message):
    __slots__ = ("read_id", "time", "channel", "action_type", "result")
    READ_ID_FIELD_NUMBER: _ClassVar[int]
    TIME_FIELD_NUMBER: _ClassVar[int]
    CHANNEL_FIELD_NUMBER: _ClassVar[int]
    ACTION_TYPE_FIELD_NUMBER: _ClassVar[int]
    RESULT_FIELD_NUMBER: _ClassVar[int]
    read_id: str
    time: float
    channel: int
    action_type: int
    result: int
    def __init__(self, read_id: _Optional[str] = ..., time: _Optional[float] = ..., channel: _Optional[int] = ..., action_type: _Optional[int] = ..., result: _Optional[int] = ...) -> None: ...

class StartRequest(_message.Message):
    __slots__ = ("acceleration_factor", "update_method", "log_interval", "stop_if_no_reads")
    ACCELERATION_FACTOR_FIELD_NUMBER: _ClassVar[int]
    UPDATE_METHOD_FIELD_NUMBER: _ClassVar[int]
    LOG_INTERVAL_FIELD_NUMBER: _ClassVar[int]
    STOP_IF_NO_READS_FIELD_NUMBER: _ClassVar[int]
    acceleration_factor: float
    update_method: str
    log_interval: int
    stop_if_no_reads: bool
    def __init__(self, acceleration_factor: _Optional[float] = ..., update_method: _Optional[str] = ..., log_interval: _Optional[int] = ..., stop_if_no_reads: bool = ...) -> None: ...

class RunMuxScanRequest(_message.Message):
    __slots__ = ("t_duration",)
    T_DURATION_FIELD_NUMBER: _ClassVar[int]
    t_duration: float
    def __init__(self, t_duration: _Optional[float] = ...) -> None: ...

class RunMuxScanResponse(_message.Message):
    __slots__ = ("nb_reads_rejected",)
    NB_READS_REJECTED_FIELD_NUMBER: _ClassVar[int]
    nb_reads_rejected: int
    def __init__(self, nb_reads_rejected: _Optional[int] = ...) -> None: ...

class BasecalledChunksRequest(_message.Message):
    __slots__ = ("batch_size", "channels")
    class Channels(_message.Message):
        __slots__ = ("value",)
        VALUE_FIELD_NUMBER: _ClassVar[int]
        value: _containers.RepeatedScalarFieldContainer[int]
        def __init__(self, value: _Optional[_Iterable[int]] = ...) -> None: ...
    BATCH_SIZE_FIELD_NUMBER: _ClassVar[int]
    CHANNELS_FIELD_NUMBER: _ClassVar[int]
    batch_size: int
    channels: BasecalledChunksRequest.Channels
    def __init__(self, batch_size: _Optional[int] = ..., channels: _Optional[_Union[BasecalledChunksRequest.Channels, _Mapping]] = ...) -> None: ...

class BasecalledReadChunkResponse(_message.Message):
    __slots__ = ("channel", "read_id", "seq", "quality_seq", "estimated_ref_len_so_far")
    CHANNEL_FIELD_NUMBER: _ClassVar[int]
    READ_ID_FIELD_NUMBER: _ClassVar[int]
    SEQ_FIELD_NUMBER: _ClassVar[int]
    QUALITY_SEQ_FIELD_NUMBER: _ClassVar[int]
    ESTIMATED_REF_LEN_SO_FAR_FIELD_NUMBER: _ClassVar[int]
    channel: int
    read_id: str
    seq: str
    quality_seq: str
    estimated_ref_len_so_far: int
    def __init__(self, channel: _Optional[int] = ..., read_id: _Optional[str] = ..., seq: _Optional[str] = ..., quality_seq: _Optional[str] = ..., estimated_ref_len_so_far: _Optional[int] = ...) -> None: ...

class BoolResponse(_message.Message):
    __slots__ = ("value",)
    VALUE_FIELD_NUMBER: _ClassVar[int]
    value: bool
    def __init__(self, value: bool = ...) -> None: ...

class DeviceInfoResponse(_message.Message):
    __slots__ = ("info", "n_channels")
    INFO_FIELD_NUMBER: _ClassVar[int]
    N_CHANNELS_FIELD_NUMBER: _ClassVar[int]
    info: str
    n_channels: int
    def __init__(self, info: _Optional[str] = ..., n_channels: _Optional[int] = ...) -> None: ...
