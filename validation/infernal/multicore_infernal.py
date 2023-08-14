# Run infernal against the input CM on all .fna files within the input directory and
# save only the tables to an output directory

import sys
from mpi4py import MPI
from subprocess import run, PIPE
import glob
import sys

in_dir = sys.argv[1]  #  input FASTA directory path
out_dir = sys.argv[2]  # where to put the output tables
cm = sys.argv[3]  # covariance model path


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# List of all fastas in the input directory
total_files = glob.glob(r"%s/*.fna" % in_dir)

# If more processes than necessary are started, exit the script
if rank >= len(total_files):
    raise SystemExit(0)

# Determine the subset of files that will be processed by one instance
count = len(total_files) // size
rem = len(total_files) % size

if rank < rem:
    start = rank * (count + 1)
    stop = start + (count + 1)
else:
    start = rank * count + rem
    stop = start + count

# Subset the main list for the current MPI instance. Each entry is the absolute path.
local_list = total_files[start:stop]

# Run infernal on every input in the subset list
for fasta in local_list:
    prod_call = run(
        [
            "cmscan",
            "--cut_ga",
            "--rfam",
            "--nohmmonly",
            "--noali",
            "--tblout",
            "%s/%s.txt" % (out_dir, fasta.split("/")[-1]),  # output file name
            "%s" % cm,
            "%s" % fasta,
        ],
        stdout=PIPE,
        stderr=PIPE,
    )

    print("%s/%s.txt" % (out_dir, fasta.split("/")[-1]))
