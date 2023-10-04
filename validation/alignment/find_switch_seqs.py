# Goes through genomes in the given dset and records riboswitch nucleotide
# sequences to files categorised by riboswitch target_name

# %%
import os
import click
import pandas as pd

from sys import argv
from glob import glob
from pathlib import Path
from Bio import SeqIO, Seq
from collections import defaultdict
from tart.utils.mpi_context import BasicMPIContext


def print(obj):
    click.echo(obj)


# %%
# Expects complete_tax_downstream.csv or
# complete_tax_downstream_mappedT_sum.csv
riboswitches_table_path = argv[1]

MAG_dir = argv[2]

# num nucleotides to capture before and after the switch
delta = int(argv[3])

dset = argv[4]

# Where to put the alignment reference fastas
out_dir = argv[5]

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
table = table[table["Dataset"] == dset]

MAG_paths = glob(f"{MAG_dir}/*.fna")
# %%
# Setup MPI to parse switch sequences
mp_con = BasicMPIContext(MAG_paths)
comm = mp_con.comm
size = mp_con.size
rank = mp_con.rank

if rank == 0:
    print("Started {} instances".format(size))

local_path_list = mp_con.generate_worker_list()

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

if mp_con.is_active:
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

seqs_arr = comm.gather(seqs_local, root=0)

# %%

# Consolidate on root thread
if rank == 0:
    print("Completed gather on 0")
    # defaultdict makes merging worker dicts trivial
    seqs_ledger = defaultdict(dict)

    for instance_dict in seqs_arr:
        for switchclass, seqs in instance_dict.items():
            seqs_ledger[switchclass].update(seqs)

    # Make the sub directory to save sequences in fasta format
    Path(f"{out_dir}/switch_seqs_delta{delta}").mkdir(parents=True, exist_ok=True)

    num_switch_classes = len(seqs_ledger)
else:
    num_switch_classes = None


# broadcast the number of distinct switch class dicts in the ledger
# so the extra threads can exit
num_switch_classes = comm.bcast(num_switch_classes, root=0)


# %%
def write_step(classname, sub_d):
    # Write riboswitch sequences to disk
    if rank > 0:
        fpath = f"{out_dir}/switch_seqs_delta{delta}/{classname}.fna"

        with open(fpath, "w") as f:
            for key, val in sub_d.items():
                f.write(">{}\n".format(key))
                f.write("{}\n".format(val))


def multithreaded_writeout(num_switch_classes, seqs_ledger):
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

    # Write riboswitch sequences to disk
    write_step(classname, sub_d)


def singlethreaded_writeout(seqs_ledger):
    if rank == 0:
        for classname, sub_d in seqs_ledger.items():
            write_step(classname, sub_d)


# %%
# %%
# Each riboswitch class sub-dictionary is sent to a worker for disk writes
# if there are enough workers for that. Otherwise the root performs the write out

if size > num_switch_classes:
    # Exit workers that are not needed
    if rank > num_switch_classes:
        raise SystemExit(0)

    if rank > 0:
        classname = None
        sub_d = None

    multithreaded_writeout(num_switch_classes, seqs_ledger)

else:
    singlethreaded_writeout(seqs_ledger)
