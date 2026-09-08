#!/usr/bin/python
"""
Some miscellaneous functions, mostly used for testing.
"""
import os
import matplotlib

#
# Headless operation.
#
# Setting STORM_ANALYSIS_HEADLESS in the environment selects a
# non-interactive matplotlib backend, which makes every pyplot.show() in
# the package a no-op. Nothing can then open a window and block, which is
# what you want when running the diagnostics, or anything else, without
# someone sitting there to close the windows.
#
# This has to happen before pyplot is imported, hence its position here.
# Saving figures is not affected, savefig() works the same under Agg.
#
if os.environ.get("STORM_ANALYSIS_HEADLESS"):
    matplotlib.use("Agg")

import matplotlib.pyplot as pyplot


class SAException(Exception):
    pass


def isHeadless():
    """
    Return True if we were asked not to open any plot windows.

    Note that this reports the request, not the outcome. matplotlib will
    also fall back to a non-interactive backend on its own when there is
    no display available.
    """
    return bool(os.environ.get("STORM_ANALYSIS_HEADLESS"))

__version__ = "2026.08.30"

# Maybe there is a builtin function that does this??
def asciiString(value):
    return str(value).encode("ascii")


def configureMatplotlib():
    """
    Configure matplotlib plots.
    """
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('legend', fontsize=10, handlelength=2)

    matplotlib.rcParams['figure.autolayout'] = True
    matplotlib.rcParams['font.size'] = 22
    matplotlib.rcParams['savefig.directory'] = os.getcwd()
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['xtick.labelsize'] = 20
    matplotlib.rcParams['xtick.major.pad'] = 10
    matplotlib.rcParams['xtick.major.size'] = 5
    matplotlib.rcParams['xtick.major.width'] = 2
    matplotlib.rcParams['xtick.top'] = 'on'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['ytick.labelsize'] = 20
    matplotlib.rcParams['ytick.major.pad'] = 10
    matplotlib.rcParams['ytick.major.size'] = 5
    matplotlib.rcParams['ytick.major.width'] = 2
    matplotlib.rcParams['ytick.right'] = 'on'    


def getData(data_path):
    import importlib.resources
    data = importlib.resources.files(__name__).joinpath(data_path)
    return str(data)


def getPath(path):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


def getPathOutputTest(fname=None):
    out_path = getPath("test/output/")
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    if fname:
        return os.path.join(out_path, fname)
    else:
        return out_path


def removeFile(fname):
    try:
        os.remove(fname)
    except OSError:
        pass

    
