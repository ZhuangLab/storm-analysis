#!/usr/bin/python
#
# Batch clustering analysis.
#
# Hazen 01/12
#

import glob
import os
import Queue
import shutil
import signal
import subprocess
import sys
import thread

# setup
if (len(sys.argv) != 3):
    print("usage: <input_directory> <channel>")
    exit()

input_directory = sys.argv[1]
channel = sys.argv[2]

src_dir = os.path.dirname(__file__)
if not (src_dir == ""):
    src_dir += "/"
    
clusters_exe = src_dir + "dbscan_analysis.py"

# find appropriate bin files
bin_files = glob.glob(input_directory + "*_alist.bin")
if(len(bin_files)==0):
    bin_files = glob.glob(input_directory + "*_list.bin")

# setup process queue
results = Queue.Queue()
def process_waiter(popen, description, que):
    try:
        popen.wait()
    finally: 
        que.put((description, popen.returncode))
process_count = 0

# start processes
procs = []
for filename in bin_files:

    # skip clustering related bin files
    if ("clusters" in filename) or ("srt" in filename):
        continue

    print("Found:", filename)

    proc = subprocess.Popen(['python', clusters_exe, filename, str(channel)])
    procs.append(proc)
    thread.start_new_thread(process_waiter, (proc, "Finished", results))
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
            proc.send_signal(signal.CTRL_C_EVENT)




