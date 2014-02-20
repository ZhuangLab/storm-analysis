#!/usr/bin/python
#
# Batch multifit analysis.
# This will start as many processes as you have files,
# so it should probably only be used on "big iron".
#
# Hazen 02/14
#

import glob
import os
import Queue
import signal
import subprocess
import sys
import thread

import sa_library.datareader as datareader

def batchAnalysis(analysis_exe, input_directory, output_directory, multi_xml, max_processes = 2):
    minimum_length = 100

    dax_files = glob.glob(input_directory + "*.dax")

    # setup process queue
    process_count = 0
    results = Queue.Queue()

    # start processes
    procs = []
    for i, file in enumerate(dax_files):

        print "Found:", file

        movie_obj = datareader.inferReader(file)
        if(movie_obj.filmSize()[2] > minimum_length):
            basename = os.path.basename(file)
            mlistname = output_directory + "/" + basename[:-4] + "_mlist.bin"
            print "  ->",mlistname

            try:
                # Wait for a process to stop before starting
                # the next one if we are at the limit.
                if(process_count >= max_processes):
                    description, rc = results.get()
                    print description
                    process_count -= 1
                proc = subprocess.Popen(['python', analysis_exe, file, mlistname, multi_xml])
                procs.append(proc)
                thread.start_new_thread(process_waiter, (proc, "Finished: " + basename, results))
                process_count += 1

            except KeyboardInterrupt:
                for proc in procs:
                    if(not proc.poll()):
                        proc.send_signal(signal.CTRL_C_EVENT)

    # wait until all the processes finish
    try:
        while(process_count>0):
            description, rc = results.get()
            print description
            process_count -= 1

    except KeyboardInterrupt:
        for proc in procs:
            if(not proc.poll()):
                proc.send_signal(signal.CTRL_C_EVENT)

def process_waiter(popen, description, que):
    try:
        popen.wait()
    finally: 
        que.put((description, popen.returncode))


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


