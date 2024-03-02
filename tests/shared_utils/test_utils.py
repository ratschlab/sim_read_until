import signal
import subprocess
import sys
import threading
import time
import pytest
from simreaduntil.shared_utils.utils import MutableValue, StoppableQueue, dill_dump, dill_load, force_eval_generator, force_eval_generator_function, get_file_content, is_sorted, num_lines_in_file, print_cmd_and_run, record_gen_fcn_waiting_time, set_signal_handler, subset_dict, tee_stdouterr_to_file, tqdm_with_name, get_some_value_from_dict
from simreaduntil.shared_utils.utils import record_gen_waiting_time

def test_is_sorted():
    import numpy as np
    for list_f in [lambda x: x, np.array]:
        assert is_sorted(list_f([7]))
        assert is_sorted(list_f([7, 9]))
        assert not is_sorted(list_f([9, 7]))
        assert is_sorted(list_f([1, 5, 8]))
        
def test_subset_dict():
    assert subset_dict({"a": 1, "b": 2}, ["a"]) == {"a": 1}
    assert subset_dict({"a": 1, "b": 2}, []) == {}
    with pytest.raises(KeyError):
        subset_dict({"a": 1, "b": 2}, ["c"])
        
def test_get_some_value_from_dict():
    assert get_some_value_from_dict({"a": 1, "b": 2}) in [1, 2] # actually 1 since order guaranteed in Python>=3.7
    
    with pytest.raises(StopIteration):
        get_some_value_from_dict({})
        
def test_num_lines_in_file(tmp_path):
    filename = tmp_path / "test.txt"
    
    def write_to_file(content):
        # overwrites file with content
        with open(filename, "w") as f:
            f.write(content)
            
    write_to_file("")
    assert num_lines_in_file(filename) == 0
    
    write_to_file("1\n2\n3")
    assert num_lines_in_file(filename) == 3
    
    # final new line, ignored
    write_to_file("1\n2\n3\n")
    assert num_lines_in_file(filename) == 3

def test_force_eval_generator():
    i = 0
    def get_numbers():
        nonlocal i
        i += 1
        yield from [1, 2, 3]
    (x for x in get_numbers())
    assert i == 0
    
    force_eval_generator(x for x in get_numbers())
    assert i == 1
    
    (x for x in force_eval_generator_function(get_numbers)())
    assert i == 2
    
def test_dill_dump_load(tmp_path):
    def get_fcn():
        i = 0
        def my_fcn():
            nonlocal i
            i += 1
            return i
        return my_fcn
    
    f = get_fcn()
    assert f() == 1
    
    save_filename = tmp_path / "fcn.pkl"
    dill_dump({"f": f, "a": 10}, save_filename)
    obj = dill_load(save_filename)
    assert obj["a"] == 10
    assert obj["f"]() == 2
    
def test_get_file_content(tmp_path):
    filename = tmp_path / "test.txt"
    with open(filename, "w") as f:
        f.write("abc\ndef")
    assert get_file_content(filename) == "abc\ndef"
    
def test_print_cmd_and_run(capfd):
    print_cmd_and_run("echo 112", verbose=True)
    assert print_cmd_and_run("echo 112", verbose=False) == "112\n"
    
    try:
        print_cmd_and_run("""
        echo 1
        >&2 echo "some error occurred"
        exit 1
        echo "never reached"
        """, verbose=False)
    except subprocess.CalledProcessError as e:
        assert "some error occurred" in e.stderr.decode()
        print("Done without errors") # previous command prints error message, to avoid confusion with real error
        
    capfd.readouterr() # clear stdout, stderr
    
def test_tqdm_with_name():
    def compute_data():
        for i in range(2):
            time.sleep(0.2)
            yield i
    for i in tqdm_with_name((i, i) for i in compute_data()):
        time.sleep(0.1)

def test_set_signal_handler():
    print_method = lambda *args, **kwargs: None
    # print_method = print

    sleep_time = 0.05
    class StoppableGen:
        def __init__(self):
            self.stop = False
            self.state = 0
            self.invalid_state = False # to signal when state is invalid, to ensure we do not exit anywhere in the code
            
        def gen_numbers(self, max_val):
            while self.state < max_val:
                self.invalid_state = False
                if self.stop: 
                    print_method("Checked stop condition, stopping")
                    break
                self.invalid_state = True
                yield self.state
                prev_state = self.state
                print_method("Sleep")
                time.sleep(sleep_time)
                print_method("after sleep")
                self.state = prev_state + 1
            self.invalid_state = False
        
    def on_keyboard_interrupt(signum, frame): 
        print_method("Keyboard interrupt, setting stop")
        gen.stop = True
        
    def raise_sigint():
        time.sleep(5*sleep_time)
        signal.raise_signal(signal.SIGINT)
    threading.Thread(target=raise_sigint).start()

    gen = StoppableGen()
    with set_signal_handler(signal.SIGINT, on_keyboard_interrupt):
        for x in gen.gen_numbers(100):
            print_method(x)
    assert not gen.invalid_state
    print_method("Continuing")
    gen.stop = False
    print_method(list(gen.gen_numbers(gen.state + 5)))
    assert not gen.invalid_state
    
def test_tee_stdouterr_to_file(tmp_path):
    with tee_stdouterr_to_file(tmp_path / "teeing", mode="w"):
        for i in range(5):
            print(f"out{i}")
            print(f"err{i}", file=sys.stderr)
            # time.sleep(0.1)
    assert (tmp_path / "teeing.out").read_text() == "".join(f"out{i}\n" for i in range(5))
    assert (tmp_path / "teeing.err").read_text() == "".join(f"err{i}\n" for i in range(5))

def test_record_gen_waiting_time():
    def produce_data():
        for i in range(5):
            time.sleep(0.1)
            yield i
        time.sleep(0.5)
        
    wait_time = MutableValue()
    for x in record_gen_waiting_time(produce_data(), wait_time):
        print(x)
        if x >= 3:
            # break early to check this is handled as well when the generator is not exhausted
            break
    # print(wait_time.value)
    assert 0.4 - 0.05 <= wait_time.value <= 0.4 + 0.05
    
def test_record_gen_fcn_waiting_time():
    def produce_data():
        for i in range(5):
            time.sleep(0.1)
            yield i
        time.sleep(0.5)
        
    def transform_data(gen):
        for x in gen:
            time.sleep(0.2)
            yield x
            
    transform_time = MutableValue()
    for x in record_gen_fcn_waiting_time(transform_data, produce_data(), transform_time):
        if x >= 3:
            # break early to check this is handled as well when the generator is not exhausted
            break
    assert 0.2*4-0.05 <= transform_time.value <= 0.2*4+0.05
    
    # as opposed to measuring the total time
    wait_time = MutableValue()
    for x in record_gen_waiting_time(transform_data(produce_data()), wait_time):
        print(x)
        if x >= 3:
            # break early to check this is handled as well when the generator is not exhausted
            break
    # print(wait_time.value)
    assert (0.2+0.1)*4-0.05 <= wait_time.value <= (0.2+0.1)*4+0.05
    
def put_item_onto_queue(queue, i, when):
    # print(f"Putting {i} onto queue at {when}")
    def f():
        queue.put(i)
        print(f"Put item '{i}' onto queue after {when}s")
    threading.Timer(when, f).start()
    
def test_StoppableQueue_NonBlockingGet():
    queue = StoppableQueue(2)
    queue.put(1)
    
    # non-blocking get
    assert queue.get(block=False) == 1
    put_item_onto_queue(queue, 2, 0.1)
    # item not yet on queue
    with pytest.raises(queue.Empty):
        queue.get(block=False)
    # block until item on queue
    assert queue.get() == 2
    
def test_StoppableQueue_FullQueueStopped():
    queue = StoppableQueue(2)
    queue.put(3)
    queue.put(4)
    # stop queue, so next get() will raise exception
    queue.stop()
    with pytest.raises(StoppableQueue.QueueStoppedException):
        queue.get()
        
def test_StoppableQueue_EmptyQueueStopped():
    queue = StoppableQueue(2)
    # empty queue that is stopped after some time
    threading.Timer(0.1, queue.stop).start()
    with pytest.raises(StoppableQueue.QueueStoppedException):
        queue.get()