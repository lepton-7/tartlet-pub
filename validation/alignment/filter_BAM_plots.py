# %%
import pickle

from pathlib import Path
from glob import glob
from mpi4py import MPI
from sys import argv
from tart.parsefuncs import plot_gen, is_interesting, bin_counts, gen_kernel

# %%
pickle_root = argv[1]
save_root = Path(argv[2])

# %%
# MPI setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# List of all pickled plots
if rank == 0:
    total_files = glob(f"{pickle_root}/**/*.p", recursive=True)

else:
    total_files = None

total_files = comm.bcast(total_files, root=0)

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

# kernel = gen_kernel(kernel_size=41, std_dev=5)

save_root = Path(save_root)
save_root.mkdir(parents=True, exist_ok=True)

bin_size = 10
for pick_file in worker_list:
    with open(pick_file, "rb") as f:
        # ([read cov, inferred frag cov, clipped cov], ends, (swtch start, end))
        cov, ends, (switch_start, switch_end) = pickle.load(f)
        readcov, infercov, clipcov = cov

    # conv_ends = np.convolve(ends, kernel, "same")

    # convalignTup = (cov, conv_ends, (switch_start, switch_end))
    alignTup = (cov, ends, (switch_start, switch_end))

    ref = Path(pick_file[:-2]).name

    core_sample = Path(pick_file).parts[-2]
    save_path = save_root.joinpath(f"{core_sample}#{ref}.png")

    if is_interesting(alignTup, windowfrac=0.2, threshtol=0.2):
        alignTup, bin_ax = bin_counts(alignTup, bin_size=bin_size)

        plot_gen(
            ref,
            alignTup,
            str(save_path),
            bin_size=bin_size,
            bin_ax=bin_ax,
        )
