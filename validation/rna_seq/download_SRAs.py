# Download sequence data using sra-tools.
# Utilises fs/scratch directory for temps
# %%
from subprocess import run, PIPE
from mpi4py import MPI
from sys import argv
from pathlib import Path

# %%
acc_file = Path(argv[1])
save_dir = argv[2]
temp_dir = argv[3]

# Make sure dirs exist
Path(save_dir).mkdir(parents=True, exist_ok=True)
Path(temp_dir).mkdir(parents=True, exist_ok=True)

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
maxTries = 3
for sra in worker_list:
    i = 0
    failed = True
    while i < maxTries and failed:
        c1 = run(
            [
                "prefetch",
                f"{sra}",
                "-O",
                f"{temp_dir}/SRAs/",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        failed = False

        if c1.returncode:
            print(f"Retrying {sra} fetch ({i})")
            print(c1.stderr)

            i += 1
            failed = True

    i = 0
    failed = True
    while i < maxTries and failed:
        c2 = run(
            [
                "fasterq-dump",
                f"{temp_dir}/SRAs/{sra}/{sra}.sra",
                "--split-files",
                "--temp",
                f"{temp_dir}",
                "--outdir",
                f"{temp_dir}",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        failed = False

        if c2.returncode:
            print(f"Retrying {sra} fastq dump ({i})")
            print(c2.stderr)

            i += 1
            failed = True

    i = 0
    failed = True
    while i < maxTries and failed:
        c3 = run(
            [
                "gzip",
                "--best",
                f"{temp_dir}/{sra}_1.fastq",
                f"{temp_dir}/{sra}_2.fastq",
            ],
            stdout=None,
            stderr=None,
        )
        failed = False

    try:  # Move the compressed archives back from temp
        movrun = run(
            [
                "mv",
                f"{temp_dir}/{sra}_1.fastq",
                f"{temp_dir}/{sra}_2.fastq",
                f"{save_dir}/",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
    except:
        print(f"mv failed for {sra}: \n{movrun.stderr}")

    try:
        delrun = run(
            [
                "rm",
                "-r",
                f"{temp_dir}/SRAs/{sra}",
            ],
            stdout=PIPE,
            stderr=PIPE,
        )

    except:
        print(f"rm failed: \n{delrun.stderr}")

    print(f"Success: {sra}")
