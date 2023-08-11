# Download sequence data using sra-tools
# %%
from subprocess import run, PIPE
from mpi4py import MPI
from sys import argv
from pathlib import Path

# %%
acc_file = Path(argv[1])
save_dir = Path(argv[2])

# Read in the list
with open(acc_file, "r") as f:
    acc_list = [x for x in f.read().splitlines()]

# %%
# MPI setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# List of SRA runs to process
total_files = acc_list

# If more processes than necessary are started, exit the script
if rank >= len(total_files):
    raise SystemExit(0)

# Determine the subset processed by one instance
count = len(total_files) // size
rem = len(total_files) % size

if rank < rem:
    start = rank * (count + 1)
    stop = start + (count + 1)
else:
    start = rank * count + rem
    stop = start + count

worker_list = total_files[start:stop]
