#!/usr/bin/python
import os


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

    
