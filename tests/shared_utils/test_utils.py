import subprocess
import time
import pytest
from simreaduntil.shared_utils.utils import dill_dump, dill_load, force_eval_generator, force_eval_generator_function, get_file_content, is_sorted, num_lines_in_file, print_cmd_and_run, subset_dict, tqdm_with_name, get_some_value_from_dict

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
