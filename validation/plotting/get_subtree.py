# Got this from: https://github.com/fredjaya/gtdb_trees

"""
Prunes a bigger tree based on a list of taxa.  

Curently used for:
    1. Subsetting group-specific taxa for treeshrink
    2. Subsetting non-empty taxa for constrained locus tree estimation
"""
# %%
import sys
from Bio import Phylo
import re
from glob import glob
from pathlib import Path

# import matplotlib # required for Phylo.draw
# %%
genpath = "../genomes"

genlist = glob(f"{genpath}/*.fna")
genlist = [".".join(Path(x).stem.split(".")[-2:]).split("_")[1] for x in genlist]

nextgen = []
for gen in genlist:
    nextgen.extend([f"GB_GCA_{gen}", f"RS_GCF_{gen}"])

# %%


def read_taxa_list(path_to_list):
    with open(path_to_list, "r") as f:
        return [line.strip() for line in f.readlines()]


def get_subtree(tree, taxa_list):
    for i, tip in enumerate(tree.get_terminals()):
        if tip.name not in taxa_list:
            try:
                tree.prune(tip)
            except ValueError:
                print(f"Failed for {i}: {tip}")
    return tree


# tree_file = sys.argv[1]
# taxa_list_file = sys.argv[2]

# tree_file = "data/bac120_r214.tree"
tree_file = "data/bac120.sp_labels.tree"
taxa_list_file = "data/taxa_list.txt"
# %%
# Read newick to prune
# tree_file = "/home/fredjaya/GitHub/gtdb_trees/data/trees/gtdb_r207_bac120_unscaled.decorated.tree"
tree = Phylo.read(tree_file, "newick")

# Read list of taxa to subset
# taxa_list_file = "/home/fredjaya/GitHub/gtdb_trees/2306_nitrospinota/nitrospinota.taxa"
taxa_list = read_taxa_list(taxa_list_file)

# %%
subtree = get_subtree(tree, taxa_list)

# Check number of subtree tips == number of taxa
# assert len([tip for tip in subtree.get_terminals()]) == len(taxa_list) // 2

# Get taxa_list basename
out = re.sub(".*/", "", taxa_list_file)
out = re.sub("\..+?$", "", taxa_list_file)
out += ".tree"
print(f"Pruned tree saved to {out}")

Phylo.write(subtree, out, "newick")
# Phylo.draw(subtree)

# %%
