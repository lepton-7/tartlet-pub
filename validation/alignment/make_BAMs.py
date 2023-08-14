# Convert riboswitch alignment SAMs to sorted BAMs

# %%
from sys import argv
from mpi4py import MPI
from subprocess import run, PIPE
from glob import glob

# %%
sam_dir = argv[1]

# %%
# MPI setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# List of all fastas in the input directory
total_files = glob(f"{sam_dir}/**/*.sam")

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

local_list = total_files[start:stop]


# %%
# SAM -> sorted BAM commands

for sam_file in local_list:
    file_base = sam_file[:-4]

    command_1_view = f"samtools view -S -b {sam_file}"
    command_2_sort = f"samtools sort {file_base}.bam -o {file_base}.sorted.bam"
    command_3_index = f"samtools index {file_base}.sorted.bam"

    with open(f"{file_base}.bam", "w") as outBAM:
        call_1 = run(command_1_view.split(" "), stdout=outBAM, stderr=PIPE)

    if call_1.returncode:
        print(call_1.stderr)

    call_2 = run(command_2_sort.split(" "), stdout=PIPE, stderr=PIPE)

    if call_2.returncode:
        print(call_2.stderr)

    call_3 = run(command_3_index.split(" "), stdout=PIPE, stderr=PIPE)

    if call_3.returncode:
        print(call_3.stderr)

    if call_1.returncode or call_2.returncode or call_3.returncode:
        print(f"Failed processing for: {sam_file}")
