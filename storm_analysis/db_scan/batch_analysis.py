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
if (len(sys.argv) != 4):
    print "usage: <input_directory> <output_directory> <channel>"
    exit()

input_directory = sys.argv[1]
output_directory = sys.argv[2]
channel = sys.argv[3]

clusters_exe = sys.path[0] + "/cl_analysis.py"

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
for file in bin_files:

    # skip clustering related bin files
    if "clusters" in file:
        continue

    print "Found:", file

    basename = os.path.basename(file)

    if 0:
        if "alist" in file:
            dir_name = output_directory + basename[:-9] + "syn"
        else:
            dir_name = output_directory + basename[:-8] + "syn"

        if os.path.exists(dir_name):
            shutil.rmtree(dir_name + "/*", ignore_errors = True)
        else:
            os.mkdir(dir_name)
        print "  ->", dir_name
    else:
        dir_name = "foo"

    proc = subprocess.Popen(['python', clusters_exe, file, dir_name + "/", str(channel)])
    procs.append(proc)
    thread.start_new_thread(process_waiter, (proc, "Finished", results))
    process_count += 1

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




