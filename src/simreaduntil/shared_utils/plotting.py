"""
Helper functions for common things in plotting
"""

import matplotlib.pyplot as plt

FIGURE_EXT = "png"
"""file extension for figures"""

def get_fig_ax(ax=None, *args, **kwargs):
    """create new figure if ax is None; otherwise, gets ax's figure"""
    if ax is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(*args, **kwargs)
    else:
        fig = ax.figure
    make_tight_layout(fig)
    return fig, ax

def make_tight_layout(fig):
    """
    Make tight-layout, compatible for different versions of matplotlib
    
    Always call this right before saving the figure because the layout may otherwise be outdated.
    """
    from matplotlib.layout_engine import TightLayoutEngine

    fig.set_layout_engine(TightLayoutEngine())
    
def ignore_tight_layout_warning():
    """
    Ignore matplotlib warning when using seaborn which calls tight_layout
    
    seaborn calls tight_layout: "UserWarning: The figure layout has changed to tight
    """
    import warnings
    warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
    
def filter_seaborn_warnings():
    """
    Filter warnings by seaborn
    """
    # import warnings
    # # python3.10/site-packages/seaborn/_oldcore.py:1119: FutureWarning: use_inf_as_na option is deprecated and will be removed in a future version. Convert inf values to NaN before operating instead.
    # # todo2: adding module="pandas" does not filter the warning, why?
    # warnings.filterwarnings("ignore", message="use_inf_as_na")
    # # python3.10/site-packages/seaborn/_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
    # warnings.filterwarnings("ignore", message="is_categorical_dtype")
    pass

def rotate_xticks(ax, rotation):
    """
    Rotate xticks of ax by rotation degrees
    
    """
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation) # only working with a warning when ticks are automatically computed, see https://stackoverflow.com/questions/43152502/how-can-i-rotate-xticklabels-in-so-the-spacing-between-each-xticklabel-is-equal
    plt.setp(ax.get_xticklabels(), ha="right", rotation=45)

def _disable_x_ticks(axis):
    """Disable x ticks on a matplotlib axis"""
    # alternative
    # axis.tick_params(axis="x", bottom=False, labelbottom=False)
    axis.set_xticks([])
    axis.set_xticklabels([])

