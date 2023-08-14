# Goes through all the derep95 MAGs and records riboswitch nucleotide
# sequences to files categorised by riboswitch target_name

# %%
from collections import defaultdict
import pandas as pd
from mpi4py import MPI
from Bio import SeqIO, Seq
from sys import argv
from glob import glob
from pathlib import Path
import os
from time import sleep

# %%
# Expects complete_tax_downstream.csv or
# complete_tax_downstream_mappedT_sum.csv
riboswitches_table_path = argv[1]

# Derep95 symlink -------------------------------------------------------------
# ENSURE THERE ARE MORE MAGs THAN RIBOSWITCH CLASSES OTHERWISE
# DOWNSTREAM MPI CODE WILL FALL APART -----------------------------------------
MAG_dir = argv[2]

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
table = table[table["Dataset"] == "Derep95"]

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
if rank >= n_records:
    raise SystemExit(0)

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

# Should have structure: {
#   "SAM" : {"SAM # 3300037155_10_Ga0395911_000841 # 6843 # 6736 # -" : "ATGC...", ...
#   },
#   "TPP" : {" TPP # ...": "GCTAGATT...", ...
#   }, ...
# }
seqs_local = defaultdict(dict)
# seqs_local = {str(rank) : rank * 50}

delta = 500  # num nucleotides to capture before and after the switch

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

        seqs_local[classname].update({rowid: switch})

# print(f"Local on rank {rank}: ", seqs_local)
# sleep(1)

seqs_arr = None
# print(f"Starting gather: {rank}")

seqs_arr = comm.gather(seqs_local, root=0)

# print(f"Completed gather: {rank}")

# %%

# Consolidate on root thread
if rank == 0:
    print("Completed gather on 0")
    # defaultdict makes merging worker dicts trivial
    seqs_ledger = defaultdict(dict)

    for instance_dict in seqs_arr:
        for switchclass, seqs in instance_dict.items():
            seqs_ledger[switchclass].update(seqs)

    # Make a directory to save sequences in fasta format
    Path("./switch_seqs_delta{}".format(delta)).mkdir(parents=True, exist_ok=True)

    num_switch_classes = len(seqs_ledger)

else:
    num_switch_classes = None


# broadcast the number of distinct switch class dicts in the ledger
# so the extra threads can exit
num_switch_classes = comm.bcast(num_switch_classes, root=0)

# %%
# Each riboswitch class sub-dictionary is sent to a worker for disk writes

# Exit workers that are not needed
if rank > num_switch_classes:
    raise SystemExit(0)

if rank > 0:
    classname = None
    sub_d = None


# Setup a dummy list to be able to subscript in the for loop and
# send each switch class to one worker thread
if not rank:
    dummy_ledger = [(key, val) for key, val in seqs_ledger.items()]

# Iterate over sub dictionaries
for idx in range(num_switch_classes):
    # Should probably refactor to use scatter instead of individual sends
    if not rank:
        key, val = dummy_ledger[idx]

        comm.send(key, dest=idx + 1, tag=10)
        comm.send(val, dest=idx + 1, tag=100)

    elif rank == idx + 1:
        classname = comm.recv(source=0, tag=10)
        sub_d = comm.recv(source=0, tag=100)

# %%
# Write riboswitch sequences to disk

if rank > 0:
    fpath = "./switch_seqs_delta{}/{}.fna".format(delta, classname)

    with open(fpath, "w") as f:
        for key, val in sub_d.items():
            f.write(">{}\n".format(key))
            f.write("{}\n".format(val))
