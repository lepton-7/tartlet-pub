# %%
import pysam
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict
from argparse import ArgumentParser
from tart.utils.mpi_context import BasicMPIContext

flagmask = 1 + 2 + 64
failmask = 3328

# CLI setup
parser = ArgumentParser(description="Extract fragment sizes from sorted BAMs.")
parser.add_argument("-i", "--input", help="Input BAMs. Can be directory or file.")
parser.add_argument("-o", "--output", help="Output .csv filepath.")


def main(parser: ArgumentParser):
    args = parser.parse_args()

    # Determine MPI context
    mp_con = BasicMPIContext()
    comm = mp_con.comm
    rank = mp_con.rank

    i_path = Path(args.input)
    ouput_path = Path(args.output)

    # Check input validity --------------------
    bam_list = []
    if i_path.is_file():
        bam_list = [f"{i_path}"]

    elif i_path.is_dir():
        bam_list = glob(f"{i_path}/*.sorted.bam")

    else:
        if not rank:
            print(f"Input neither a file nor a directory: {i_path}")

        SystemExit(0)

    mp_con.set_full_list(bam_list)
    worker_list = mp_con.generate_worker_list()

    # bam_list = [
    #     "/fs/scratch/PDS0325/rna_seq/e_coli_alignment_20231113/SRR7154636.sorted.bam",
    #     "/fs/scratch/PDS0325/rna_seq/e_coli_alignment_20231113/SRR7154637.sorted.bam",
    # ]

    df_recs = []
    for bampath in worker_list:
        sra = Path(bampath).stem.split(".")[0]

        print(f"Processing {sra} on worker {rank}.")

        len_dict = defaultdict(int)
        bam = pysam.AlignmentFile(bampath, "rb")

        # List of contigs in the bam
        ref_list = [
            idxstats.contig
            for idxstats in bam.get_index_statistics()
            if idxstats.total > 0
        ]

        # Go through each contig and make sure to fully traverse the BAM
        for ref in ref_list:
            for pys_read in bam.fetch(contig=ref, multiple_iterators=True):
                if (pys_read.flag & flagmask == flagmask) and (
                    pys_read.flag & failmask == 0
                ):
                    len_dict[str(abs(pys_read.template_length))] += 1

        df = pd.DataFrame.from_records(list(len_dict.items()))
        df.columns = ["size", sra]

        df_recs.append(df)

    # DF merge on worker
    df = pd.DataFrame({"size": []})
    for temp in df_recs:
        df = pd.merge(df, temp, how="outer", on="size")

    df_arr = comm.gather(df, root=0)

    if rank == 0:
        if df_arr is None:
            raise ValueError("Gather failed on root.")
        print("Completed gather.")

        print("Starting merge.")
        final_df = pd.DataFrame({"size": []})
        for df_inst in df_arr:
            final_df = pd.merge(final_df, df_inst, how="outer", on="size")
        print("Finished merge.")

        final_df.fillna(0, inplace=True)
        print("Sorting.")
        final_df.sort_values("size", inplace=True, key=pd.to_numeric)

        print("Completing output.")
        final_df.to_csv("sizes.csv", index=False)

    SystemExit(0)


if __name__ == "__main__":
    main(parser)
