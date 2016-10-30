#!/usr/bin/python
#
# Handle loading the correct C library.
#
# Hazen 11/14
#

import ctypes
import sys
import os
import re

import storm_analysis

def loadCLibrary(package, library):

    # Something like /usr/lib/python3.5/site-packages/storm_analysis
    module_path = os.path.dirname(os.path.dirname(os.path.abspath(storm_analysis.__file__)))

    # Something like /usr/lib/python3.5/site-packages/storm_analysis/sa_library/
    lib_path = os.path.join(module_path, package.replace(".", os.path.sep))
    files = os.listdir(lib_path)

    lib_extension = "so"
    if (sys.platform == "win32"):
        lib_extension = "dll"
                
    # Something like '_ia_utilities.*\.so'
    r = re.compile('{}.*\.{}'.format(library, lib_extension))

    lib_filename = list(filter(r.match, files))
    
    if len(lib_filename) < 1:
        raise Exception("Can't find the library {} in the module {} "
                        "located in the storm_analysis package at {}".format(library, package, module_path))
    
    # Something like _ia_utilities.cpython-35m-x86_64-linux-gnu.so
    lib_filename = lib_filename[0]
    
    return ctypes.cdll.LoadLibrary(os.path.join(lib_path, lib_filename))

#
# The MIT License
#
# Copyright (c) 2014 Zhuang Lab, Harvard University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
