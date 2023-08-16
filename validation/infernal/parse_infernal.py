# %%
import pandas as pd

# %%
tax_cache = {}

# %%
# infernal output tables and dataset identifiers
files = {
    "ecoli": (r"./results/Ecoli.NC_000913.3.fna.txt", "e_coli"),
    "bsub": (r"./results/Bsubtilis168.GCF_000009045.1.fna.txt", "b_sub_168"),
    "pfluor": (r"./results/Pfluorescens.GCA_900215245.1.fna.txt", "p_fluor"),
    "paeru": (r"./results/Paeruginosa.GCA_000006765.1.fna.txt", "p_aeru"),
    "mtuber": (r"./results/MtuberculosisH37Rv.GCA_000195955.2.fna.txt", "m_tuber"),
    "selon": (r"./results/Selongatus.GCA_000012525.1.fna.txt", "s_elon"),
}

table_dir = "/users/PDS0325/sachitk26/packages/tart/validation/tables"


def parse_infernal(process_file, name):
    rows = []
    # This is so that extra data that is not in the infernal table can be
    # added to the row data and parsed at the same time
    additionalCols = r" %s %s %s %s %s %s %s %s %s %s %s"

    with open(process_file, "r") as f:
        MAG_acc = " "

        # infernal puts the query file (which has its MAG accession ID) at
        # the end of the riboswitch tables, so the easy way to associate
        # MAG IDs with their riboswitches is to just reverse the results
        # and append it to the results until the next table and a new MAG ID is detected
        for line in reversed(f.readlines()):
            li = line.strip()

            # ignore everything except the results rows and add taxonomy data
            if not li.startswith("#"):
                # format the taxonomy string to determine class
                try:
                    tax_entry = tax_cache[MAG_acc]
                    # Get rid of the species name as the space interfers with downstream information
                    # consolidation
                    splits = tax_entry[0].split(";")
                    if splits[-1].split("__")[0] == "s":
                        tax_entry[0] = ";".join(splits[:-1])

                except KeyError:
                    tax_entry = ["NA", "NA", "NA"]

                    # taxonomy string looks like
                    # d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;f__UBA4665;g__UBA5150;s__
                    # first split by ";" and select element 2: c__Saccharimonadia
                    # then split by "__" and select element 1: Saccharimonadia

                try:
                    MAG_dom = tax_entry[0].split(";")[0].split("__")[1]
                except IndexError:
                    MAG_dom = "NA"

                try:
                    MAG_phylum = tax_entry[0].split(";")[1].split("__")[1]
                except IndexError:
                    MAG_phylum = "NA"

                try:
                    MAG_class = tax_entry[0].split(";")[2].split("__")[1]
                except IndexError:
                    MAG_class = "NA"

                try:
                    MAG_order = tax_entry[0].split(";")[3].split("__")[1]
                except IndexError:
                    MAG_order = "NA"

                try:
                    MAG_fam = tax_entry[0].split(";")[4].split("__")[1]
                except IndexError:
                    MAG_fam = "NA"

                try:
                    MAG_genus = tax_entry[0].split(";")[5].split("__")[1]
                except IndexError:
                    MAG_genus = "NA"

                rows.append(
                    line.rstrip()
                    + additionalCols
                    % (
                        MAG_acc,
                        name,
                        MAG_dom,
                        MAG_phylum,
                        MAG_class,
                        MAG_order,
                        MAG_fam,
                        MAG_genus,
                        *tax_entry,
                    )  # tax_entry has 3 elements
                )

            # start of results from a new query MAG
            if li.startswith("# Query file:"):
                # Extract the MAG ID from the big query file path. Example:
                # Query file: /users/PDS0325/sachitk26/Emerge_MAGs_v1/Derep95/20110600_P1M_23.fna
                # first at every /, then take the last bit: 20110600_P1M_23.fna
                # and slice it to remove the ".fna".
                MAG_acc = line.rstrip().split("/")[-1][:-4]

    table = []

    for row in rows:
        # Taking the strings from the infernal output and putting them into
        # a table to turn into a pd df.
        row_contents = row.split()

        # Because data in the "description of target" column has spaces,
        # the split() commands puts each word as its own element
        # into row_contents. To pull the relevant info out of the array,
        # slice and include in a way that doesnt touch
        # indices occupied by the description.
        table.append(
            [*row_contents[:16], *row_contents[-len(additionalCols.split(" ")) + 1 :]]
        )

    # Took most of this from the infernal results table header to pass onto the df
    col_lab = [
        "target_name",
        "accession",
        "query_name",
        "accession",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "E-value",
        "MAG_accession",
        "Dataset",
        "MAG_domain",
        "MAG_phylum",
        "MAG_class",
        "MAG_order",
        "MAG_family",
        "MAG_genus",
        "full_taxonomy",
        "completeness",
        "contamination",
    ]

    df = pd.DataFrame(
        table,
        columns=col_lab,
    )

    # I dont think these columns are informative. But then again I also have 0 clue what they mean
    df = df.drop(columns=["mdl_from", "mdl_to", "mdl"], axis="columns")

    return df


# %%

dfs = []
for dataset in files:
    path, name = files[dataset]
    dfs.append(parse_infernal(path, name))

# %%

filename = "inf_results"
cumulative = pd.concat(dfs, ignore_index=True)
cumulative.to_csv(
    f"{table_dir}/{filename}.csv",
    index=False,
)

# %%
