# Goes through derep95 MAGs, finds riboswitch sequences with a buffer around the sequence
# and returns a dict for indices of the full switch + buffer sequences in the contig

# %%
from collections import defaultdict
import pandas as pd
from mpi4py import MPI
from Bio import SeqIO, Seq
from sys import argv
from glob import glob
import os
import json

# %%
# Expects complete_tax_downstream.csv or
# complete_tax_downstream_mappedT_sum.csv
riboswitches_table_path = argv[1]

# Derep95 symlink -------------------------------------------------------------
# ENSURE THERE ARE MORE MAGs THAN RIBOSWITCH CLASSES OTHERWISE
# DOWNSTREAM MPI CODE WILL FALL APART -----------------------------------------
MAG_dir = argv[2]

# num nucleotides to capture before and after the switch
delta = int(argv[3])

dset = argv[4]

# Read in the results table
table = pd.read_csv(riboswitches_table_path)

# Drop all but the relevant columns
table = table[
    [
        "target_name",
        "query_name",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "MAG_accession",
        "Dataset",
    ]
]
# table = table[table["Dataset"] == dset]

MAG_paths = glob("{}/*.fna".format(MAG_dir))
# %%
# Setup MPI to parse switch sequences

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    print("Started {} instances".format(size))

# number of derep95 MAGs
n_records = len(MAG_paths)


# If more processes than necessary are started, exit the script
skip = False
if rank >= n_records:
    skip = True
    # raise SystemExit(0)

count = n_records // size
rem = n_records % size

# Determine the subset of MAGs that will be processed by one instance
if rank < rem:
    startMPI = rank * (count + 1)
    stopMPI = startMPI + (count + 1)
else:
    startMPI = rank * count + rem
    stopMPI = startMPI + count


local_path_list = MAG_paths[startMPI:stopMPI]

# %%
# Setup the local dictionary that stores the sequences
seqs_local = defaultdict()
if not skip:
    for MAG_path in local_path_list:  # iterate over derep95 MAGs
        MAGDict = {x.id: str(x.seq) for x in SeqIO.parse(MAG_path, "fasta")}
        subset = table[table["MAG_accession"] == os.path.split(MAG_path)[-1][:-4]]

        for _, row in subset.iterrows():  # iterate over MAG riboswitches
            # Infernal start and stop entries are relative to canonical 5' -> 3';
            # need to account for that when slicing the sequence
            if row["strand"] == "+":
                start = int(row["seq_from"])
                end = int(row["seq_to"])

            elif row["strand"] == "-":
                start = int(row["seq_to"])
                end = int(row["seq_from"])

            else:
                print("Yikes")  # should not be possible

            contigseq = MAGDict[row["query_name"]]

            # Adjust bounds with delta
            start -= delta
            end += delta

            # Validate bounds
            if start < 0:
                start = 0
            if end > len(contigseq):
                end = len(contigseq)

            # Find query seq
            switch = contigseq[start:end]

            # switch is on the (-) strand
            if row["strand"] == "-":
                switch = Seq.reverse_complement(switch)

            # Add the sequence to the dict
            classname = row["target_name"]
            qname = row["query_name"]
            frm = row["seq_from"]
            to = row["seq_to"]
            strand = row["strand"]

            # No spaces to ensure the entire string is recognised as the ID
            rowid = f"{classname}#{qname}#{frm}#{to}#{strand}"

            seqs_local.update({rowid: (start, end)})

        print(seqs_local)

seqs_arr = None

seqs_arr = comm.gather(seqs_local, root=0)

# print(f"Completed gather: {rank}")

# %%

# Consolidate on root thread
if rank == 0:
    print("Completed gather on 0")
    # defaultdict makes merging worker dicts trivial
    seqs_ledger = defaultdict()

    for instance_dict in seqs_arr:
        seqs_ledger.update(instance_dict)

    with open("rowid_to_bounds.json", "w") as f:
        json.dump(seqs_ledger, f)

else:
    raise SystemExit(0)
