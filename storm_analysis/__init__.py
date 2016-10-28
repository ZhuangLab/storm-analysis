#!/usr/bin/python
import os


def get_data(data_path):
    import pkg_resources
    data = pkg_resources.resource_filename(__name__, data_path)
    return data


def get_path(path):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


def get_path_output_test(fname=None):
    out_path = get_path("test/output/")
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    if fname:
        return os.path.join(out_path, fname)
    else:
        return out_path