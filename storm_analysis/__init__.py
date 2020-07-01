#!/usr/bin/python
"""
Some miscellaneous functions, mostly used for testing.
"""
import os
import matplotlib
import matplotlib.pyplot as pyplot


class SAException(Exception):
    pass

__version__ = "2020.07.01"

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
    matplotlib.rcParams['savefig.directory'] = os.chdir(os.getcwd())
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
    import pkg_resources
    data = pkg_resources.resource_filename(__name__, data_path)
    return data


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

    
