# %%
import pysam
import pandas as pd
import matplotlib.pyplot as plt

from glob import glob
from pathlib import Path
from collections import defaultdict
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
    "/fs/scratch/PDS0325/rna_seq/e_coli_alignment_20231113/SRR7154636.sorted.bam",
    "/fs/scratch/PDS0325/rna_seq/e_coli_alignment_20231113/SRR7154637.sorted.bam",
]

df_recs = []
for bampath in bam_list:
    sra = Path(bampath).stem.split(".")[0]

    len_dict = defaultdict(int)
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
                len_dict[str(abs(pys_read.template_length))] += 1
                i += 1

                if i % 200_000 == 0:
                    print(f"Processed {i} reads")

    df = pd.DataFrame.from_records(list(len_dict.items()))
    df.columns = ["size", sra]

    df_recs.append(df)


# %%
df = pd.DataFrame({"size": []})

for temp in df_recs:
    df = pd.merge(df, temp, how="outer", on="size")

df.fillna(0, inplace=True)
df.sort_values("size", inplace=True, key=pd.to_numeric)

# %%
df.to_csv("sizes.csv", index=False)

# %%
