#!/usr/bin/env python
"""
Batch clustering analysis.

Hazen 01/12
"""

import glob
import os
import multiprocessing
import shutil
import signal
import subprocess
import sys
import threading


def batchAnalysis(input_directory, channel, eps = 40, mc = 10, min_size = 50):

    src_dir = os.path.dirname(__file__)
    if not (src_dir == ""):
        src_dir += "/"
    
    clusters_exe = src_dir + "dbscan_analysis.py"

    # find appropriate bin files
    bin_files = glob.glob(input_directory + "*_alist.bin")
    if(len(bin_files)==0):
        bin_files = glob.glob(input_directory + "*_list.bin")

    # setup process queue
    results = multiprocessing.Queue()
    def process_waiter(popen, description, que):
        try:
            popen.wait()
        finally: 
            que.put((description, popen.returncode))

    # start processes
    process_count = 0
    procs = []
    for filename in bin_files:

        # skip clustering related bin files
        if ("clusters" in filename) or ("srt" in filename):
            continue

        print("Found:", filename)

        proc = subprocess.Popen(['python', clusters_exe,
                                 "--bin", filename,
                                 "--channel", str(channel),
                                 "--eps", str(eps),
                                 "--mc", str(mc),
                                 "--min_size", str(min_size)])
        procs.append(proc)
        t = threading.Thread(target = process_waiter,
                             args = (proc, "Finished", results))
        t.daemon = True
        t.start()
        process_count += 1

    # wait until all the processes finish
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


if (__name__ == "__main__"):

    import argparse

    parser = argparse.ArgumentParser(description = 'Batch DBSCAN clustering.')

    parser.add_argument('--dir', dest='directory', type=str, required=True,
                        help = "The name of the directory containing the Insight3 format files to analyze.")
    parser.add_argument('--channel', dest='channel', type=int, required=True,
                        help = "Which channel (or category) to use for clustering.")
    parser.add_argument('--eps', dest='epsilon', type=float, required=False, default=40,
                        help = "The DBSCAN epsilon parameters in nanometers. The default is 40nm.")
    parser.add_argument('--mc', dest='mc', type=int, required=False, default=10,
                        help = "The DBSCAN mc parameter. The default is 10.")
    parser.add_argument('--min_size', dest='min_size', type=int, required=False, default=50,
                        help = "The minimum cluster size to include when calculating cluster statistics. The default is 50.")

    args = parser.parse_args()

    batchAnalysis(args.directory, args.channel, args.epsilon, args.mc, args.min_size)

