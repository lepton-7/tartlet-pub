# %%
import pysam
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
from argparse import ArgumentParser
from tart.utils.mpi_context import BasicMPIContext

# %%

flagmask = 1 + 2 + 64
failmask = 3328

# %%


# def main():
# Determine MPI context
# mp_con = BasicMPIContext()
# comm = mp_con.comm
# rank = mp_con.rank

# bam_list = []
bam_list = [
    "/fs/scratch/PDS0325/rna_seq/e_coli_alignment_20231113/SRR7154636.sorted.bam"
]
len_list = []
for bampath in bam_list:
    len_list = []
    bam = pysam.AlignmentFile(bampath, "rb")

    # List of contigs in the bam
    ref_list = [
        idxstats.contig for idxstats in bam.get_index_statistics() if idxstats.total > 0
    ]

    i = 0
    for ref in ref_list:
        for pys_read in bam.fetch(contig=ref, multiple_iterators=True):
            if (pys_read.flag & flagmask == flagmask) and (
                pys_read.flag & failmask == 0
            ):
                len_list.append(abs(pys_read.template_length))
                i += 1

                if i % 100000 == 0:
                    print(f"Processed {i} reads")


# %%
df = pd.DataFrame({"sizes": len_list})
df.to_csv("sizes.csv", index=False)
