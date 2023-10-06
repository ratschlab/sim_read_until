from matplotlib import pyplot as plt
import pytest

from simreaduntil.shared_utils.plotting import make_tight_layout
from simreaduntil.simulator.utils import format_percentage, in_interval, new_thread_name


def test_new_thread_name():
    name1 = new_thread_name(template_str="ontsim-{}")
    name2 = new_thread_name(template_str="ontsim-{}")
    assert name1.startswith("ontsim-")
    assert name2.startswith("ontsim-")
    assert name1 != name2
    
def test_in_interval():
    assert in_interval(0.5, (0, 1))
    assert not in_interval(0.5, (0, 0.4))
    
def test_format_percentage():
    assert format_percentage(0.51034) == "51.03"
    
def test_make_tight_layout():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    make_tight_layout(fig)