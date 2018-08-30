#!/usr/bin/env python
"""
Handle loading the correct C library.

Hazen 11/14
"""

import ctypes
import sys
import os
import re

import storm_analysis

def loadCLibrary(library_filename):

    #
    # This assumes that all the code that will call C modules is one level
    # below the root storm_analysis directory.
    #
    # c_lib_path is something like:
    #    /usr/lib/python3.5/site-packages/storm_analysis
    #
    c_lib_path = os.path.dirname(os.path.abspath(storm_analysis.__file__))

    # All the C libraries are in the c_libraries directory.
    c_lib_path = os.path.join(c_lib_path, "c_libraries")

    # Windows.
    if (sys.platform == "win32"):

        library_filename += '.dll'

        # Try to load the library without fiddling with the Windows DLL search 
        # first. Do this so that if the user wanted to use their own versions 
        # of FFTW or LAPACK these will get loaded instead of the ones in 
        # storm-analysis.
        try:
            # This suppresses the Windows DLL missing dialog.
            ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
            return ctypes.cdll.LoadLibrary(os.path.join(c_lib_path, library_filename))
        except WindowsError:
            # Unsuppress the Windows DLL missing dialog.
            ctypes.windll.kernel32.SetErrorMode(0)

            # Push the storm-analysis directory into the DLL search path.
            if (sys.version_info > (3, 0)):
                ctypes.windll.kernel32.SetDllDirectoryW(c_lib_path)
            else:
                ctypes.windll.kernel32.SetDllDirectoryW(unicode(c_lib_path))

            # Try to load the library.
            c_lib = ctypes.cdll.LoadLibrary(os.path.join(c_lib_path, library_filename))

            # Restore the Windows DLL search path.
            ctypes.windll.kernel32.SetDllDirectoryW(None)

            return c_lib

    # OS-X.
    elif (sys.platform == "darwin"):
        library_filename = 'lib' + library_filename
        library_filename += '.dylib'
        return ctypes.cdll.LoadLibrary(os.path.join(c_lib_path, library_filename))
        
    # Linux.
    else:
        library_filename = 'lib' + library_filename
        library_filename += '.so'
        return ctypes.cdll.LoadLibrary(os.path.join(c_lib_path, library_filename))


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
