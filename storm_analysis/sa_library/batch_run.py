#!/usr/bin/env python
"""
Run multiple python scripts in parallel.

Hazen 02/17
"""

import os
import multiprocessing
import signal
import subprocess
import sys
import threading


def batchRun(cmd_lines, max_processes = 2):

    # Setup process queue.
    process_count = 0
    results = multiprocessing.Queue()

    # Start processes.
    procs = []
    for cmd_line in cmd_lines:

        try:
            # Wait for a process to stop before starting
            # the next one if we are at the limit.
            if(process_count >= max_processes):
                description, rc = results.get()
                print(description)
                process_count -= 1
            proc = subprocess.Popen(cmd_line,
                                    env = os.environ.copy())
            procs.append(proc)
            t = threading.Thread(target = processWaiter,
                                 args = (proc, "Finished '" + " ".join(cmd_line) + "'", results))
            t.daemon = True
            t.start()
            process_count += 1

        except KeyboardInterrupt:
            for proc in procs:
                if(not proc.poll()):
                    if (sys.platform == "win32"):
                        proc.send_signal(signal.CTRL_C_EVENT)
                    else:
                        proc.send_signal(signal.SIGINT)
            break

    # Wait until all the processes finish.
    try:
        while(process_count>0):
            description, rc = results.get()
            print(description)
            process_count -= 1

    except KeyboardInterrupt:
        for proc in procs:
            if(not proc.poll()):
                if (sys.platform == "win32"):
                    proc.send_signal(signal.CTRL_C_EVENT)
                else:
                    proc.send_signal(signal.SIGINT)
                    

def processWaiter(popen, description, que):
    try:
        popen.wait()
    finally: 
        que.put((description, popen.returncode))


#
# The MIT License
#
# Copyright (c) 2017 Zhuang Lab, Harvard University
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


