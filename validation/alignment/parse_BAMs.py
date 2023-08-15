# %%
import pysam
import json, pickle
from glob import glob

from mpi4py import MPI
from sys import argv
from pathlib import Path
from tart.parsefuncs import generate_plot_data, plot_gen

# %%
bam_dir = argv[1]

save_root = Path(argv[2])
save_root.mkdir(parents=True, exist_ok=True)

# Set these to determine outputs ----------------------------------------------
minfragcov = 8

outPlots = False
outPickles = True
# -----------------------------------------------------------------------------

# Import the dict that stores MAG contig positions of the alignment references
with open("./rowid_to_bounds.json", "r") as f:
    refBounds = json.load(f)

# %%
# MPI setup
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# List of all sorted BAMS to process
total_files = glob(f"{bam_dir}/**/*.sorted.bam")

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
for bam_path in worker_list:
    bam = pysam.AlignmentFile(bam_path, "rb")

    # This looks something like this:
    # /users/PDS0325/sachitk26/ribo_acclim/metatranscriptomics/
    # switch_alignment/switch_seqs_delta500/alignments_full/AdoCbl_riboswitch/AdoCbl_riboswitch.MainAutochamber.201707_E_2_20to24.sorted.bam
    bam_path = Path(bam_path)
    bam_split = bam_path.parts

    # Looks like:
    # save_dir/AdoCbl_riboswitch/AdoCbl_riboswitch.MainAutochamber.201707_E_2_20to24
    #
    # The [:-11] removes the .sorted.bam suffix from the path, but Path objects are
    # not subscriptable, so it needs to be typecasted to str, sliced, then converted
    # back to a Path. Don't @ me this works just fine.
    save_dir = Path(str(save_root.joinpath(*bam_split[-2:]))[:-11])

    # Make sure the directory exists; create if not.
    save_dir.mkdir(parents=True, exist_ok=True)

    outDict = generate_plot_data(bam, refBounds)

    for ref, alignTup in outDict.items():
        start, end = alignTup[2][0], alignTup[2][1]
        readcov = alignTup[0][0]

        if outPlots and max(readcov[start:end]) >= minfragcov:
            save_path = save_dir.joinpath(f"{ref}.png")
            plot_gen(ref, alignTup, save_path, buff=40)

        if outPickles and max(readcov[start:end]) >= minfragcov:
            save_path = save_dir.joinpath(f"{ref}.p")
            with open(save_path, "wb") as f:
                pickle.dump(alignTup, f)
