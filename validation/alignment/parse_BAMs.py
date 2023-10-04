# %%
import pysam
import json, pickle

from glob import glob
from sys import argv
from pathlib import Path
from tart.utils.mpi_context import BasicMPIContext
from tart.utils.parsefuncs import generate_plot_data, plot_gen

# %%
bam_dir = argv[1]

save_root = Path(argv[2])
save_root.mkdir(parents=True, exist_ok=True)

bounds_dir = argv[3]

# Set these to determine outputs ----------------------------------------------
minfragcov = 8

outPlots = False
outPickles = True
# -----------------------------------------------------------------------------

# Import the dict that stores MAG contig positions of the alignment references
with open(f"{bounds_dir}/rowid_to_bounds.json", "r") as f:
    refBounds = json.load(f)

# %%
# List of all sorted BAMS to process
total_files = glob(f"{bam_dir}/**/*.sorted.bam")

# MPI setup
mp_con = BasicMPIContext(total_files)
worker_list = mp_con.generate_worker_list()

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
