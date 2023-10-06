"""
Utility functions to merge plt axes into one axis includings lines (ax.plot), bars (ax.bar, ax.hist)
"""

from pathlib import Path
from matplotlib import pyplot as plt
from simreaduntil.shared_utils.plotting import make_tight_layout
from simreaduntil.shared_utils.utils import dill_dump, dill_load

def copy_axis_props(ax_from, ax_to, **kwargs):
    """
    Copy the line from ax_from to ax_to
    """
    ax_to.set_title(ax_from.get_title())
    ax_to.set_xlabel(ax_from.get_xlabel())
    ax_to.set_xscale(ax_from.get_xscale())
    ax_to.set_ylabel(ax_from.get_ylabel())
    ax_to.set_yscale(ax_from.get_yscale())
    
def compare_axes(*axes, ignore_empty=True):
    """Compares properties of axes to see if they are different"""
    res = {
        "title": [ax.get_title() for ax in axes],
        "xlabel": [ax.get_xlabel() for ax in axes],
        "xscale": [ax.get_xscale() for ax in axes],
        "ylabel": [ax.get_ylabel() for ax in axes],
        "yscale": [ax.get_yscale() for ax in axes],
    }
    new_res = {}
    for key, val in res.items():
        if ignore_empty:
            val = [x for x in val if x != ""]
        if len(set(val)) > 1:
            new_res[key] = set(val)
    return new_res

def copy_line(line, ax, **kwargs):
    """
    Copy the line to the axis
    
    A matplotlib artist (e.g. line) can only belong to one axis at a time, so we have to either remove
    it if we want to put it onto another axis, or copy the data.
    
    Use default kwargs by setting them to None, e.g. label=None
    """
    # ax.add_line(Line2D(*line.get_data(), label=line.get_label(), **kwargs))
    # better since the axis autoscales
    extra_args = {
        "label": line.get_label(),
        # "color": line.get_color(), # let it cycle through colors itself
        "linestyle": line.get_linestyle(),
        "linewidth": line.get_linewidth(),
        "marker": line.get_marker(),
        "markersize": line.get_markersize(),
        "markerfacecolor": line.get_markerfacecolor(),
        "markeredgecolor": line.get_markeredgecolor(),
        "markeredgewidth": line.get_markeredgewidth(),
        "markevery": line.get_markevery(),
        "alpha": line.get_alpha(),
        "zorder": line.get_zorder(),
    }
    # alternatively can use Line2D.update_from
    extra_args.update(kwargs)
    ax.plot(*line.get_data(), **extra_args)
    
# note: copying a histogram is not too useful, since it is a collection of patches and does not allow rebinning
# when plotting side-by-side with another dataset

def save_fig_and_pickle(fig, filename):
    """
    Save the figure and pickle it for later loading, so it can be merged with another figure
    
    Uses dill
    """
    fig.savefig(filename)
    filename = Path(filename).with_suffix('.dill')
    
    # put into subdirectory of where figure pdf/png is saved
    filename = filename.parent / "pickled_figures" / filename.name
    filename.parent.mkdir(exist_ok=True)
    
    dill_dump(fig, filename)
    
if __name__ == "__main__":
    import numpy as np

    # Create some example data for multiple plots
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x)
    y2 = np.cos(x)
    y3 = np.tan(x)

    fig0, (ax1, ax2) = plt.subplots(ncols=2)
    ax1.plot(x, y1, label="sin", color="g", linestyle="--")
    ax1.set_title("figure title")
    # ax1.set_xscale("log")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim([0, 15])
    ax1.legend()
    ax2.plot(x, y2, label="cos")
    ax2.set_xlim([0, 20])
    ax2.set_xlabel("x2")
    # plt.close(fig0)

    compare_axes(ax1, ax2)
    compare_axes(ax1, ax2, ignore_empty=False)

    fig, ax = plt.subplots()
    copy_axis_props(ax1, ax)
    line1 = ax1.get_lines()[0]
    line2 = ax2.get_lines()[0]
    copy_line(line1, ax)
    copy_line(line2, ax, color="r")
    ax.legend()
    ax.autoscale()
    make_tight_layout(ax.figure)
    
    
    # save and load
    
    fig1, ax1 = plt.subplots()
    ax1.plot(x, y1, label="sin", color="g", linestyle="--")
    ax1.set_title("figure title")
    # ax1.set_xscale("log")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim([0, 15])
    ax1.legend()
    fig2, ax2 = plt.subplots()
    ax2.plot(x, y2, label="cos")
    ax2.set_xlim([0, 20])
    ax2.set_xlabel("x2")

    figure_filename1 = Path('figure1.png')
    figure_filename2 = Path('figure2.png')
    save_fig_and_pickle(fig1, figure_filename1)
    save_fig_and_pickle(fig2, figure_filename2)
    del fig1, fig2, ax1, ax2
    fig1 = dill_load(figure_filename1.with_suffix('.dill'))
    fig2 = dill_load(figure_filename2.with_suffix('.dill'))
    ax1 = fig1.axes[0]
    ax2 = fig2.axes[0]

    fig, ax = plt.subplots()
    copy_axis_props(ax1, ax)
    line1 = ax1.get_lines()[0]
    line2 = ax2.get_lines()[0]
    copy_line(line1, ax)
    copy_line(line2, ax, color="r")
    ax.legend()
    ax.autoscale()
    make_tight_layout(ax.figure)
