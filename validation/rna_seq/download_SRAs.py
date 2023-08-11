# Download sequence data using sra-tools
# %%
from subprocess import run, PIPE, CalledProcessError
from mpi4py import MPI
from sys import argv
from pathlib import Path

# %%
acc_file = Path(argv[1])
save_dir = argv[2]

# save_dir = "/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/e_coli"
# acc_file = "/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/e_coli_sastry2019_acc.txt"

# sra = "SRR8164486"

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

# %%
for sra in worker_list:
    retries = 0
    failed = True
    while retries < 5 and failed:
        try:
            c1 = run(
                [
                    "prefetch",
                    f"{sra}",
                    "-O",
                    "SRAs/"
                ],
                stdout=PIPE,
                stderr=PIPE,
                check=True)

            c2 = run(
                [
                    "fasterq-dump",
                    f"SRAs/{sra}/{sra}.sra",
                    "--split-files",
                    "--outdir",
                    f"{save_dir}"
                ],
                stdout=PIPE,
                stderr=PIPE,
                check=True
            )
            failed = False
            print(f"Success: {sra} ")

        except CalledProcessError as e:
            retries += 1
            failed = True
            print(f"Failed to process {sra}. Retrying {retries+1}...")
            print(f"{sra}: {e.stderr}")

            if retries == 5:
                print(f"Skipping {sra}")

