# Collects data for the big results figure from each peak plot cluster stats
#  file from the validation set

# %%
from typing import Literal
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict


# %%
def classname(rowid: str):
    return rowid.split("#")[0]


mean_t = 0.05
var_t = 0.05

# %%

# -----------------------------------------------
# Exclude these datasets
# -----------------------------------------------
excluded_dsets = [
    "a_baum",
    "c_basil",
    "p_fluor",
    "p_aeru",
    "p_salmonis",
]

# Find all cluster stats files
cluster_paths = [
    Path(x)
    for x in glob("../alignment/outputs/*/plots/cluster_stats.csv")
    if Path(x).parent.parent.stem not in excluded_dsets
]

peak_paths = [
    Path(x)
    for x in glob("../alignment/outputs/*/plots/peak_log.csv")
    if Path(x).parent.parent.stem not in excluded_dsets
]

# Ingest the table with literature evidence.
lit_table_path = "./data/big_fig/lit_findings.csv"
lit_table = pd.read_csv(lit_table_path)

switch_table_path = "../tables/inf_results.csv"


# %%
def make_rowid(r: pd.Series):
    return f'{r["target_name"]}#{r["query_name"]}#{r["seq_from"]}#{r["seq_to"]}#{r["strand"]}'


def lit_agreement(rowid: str, tartres: int):
    subdf = lit_table[lit_table["rowid"] == rowid]

    if len(subdf) == 0:
        return None

    if tartres and subdf["lit_active"].iloc[0]:
        return "active_agree"

    elif tartres and not subdf["lit_active"].iloc[0]:
        return "active_disagree"

    elif not tartres and subdf["lit_active"].iloc[0]:
        return "incon_disagree"

    elif not tartres and subdf["lit_incon"].iloc[0]:
        return "incon_agree"


# %%
# Go through each peak_log and find how many transcriptomes contribute to peaks for each dataset

# TODO: Talk to Sarah about whether it's a better to log number of transcriptomes input instead.
sra_dict = {}

for ppath in peak_paths:

    dset = ppath.parent.parent.stem

    df = pd.read_csv(ppath)
    sra_dict[dset] = len(pd.unique(df["transcriptome"]))


# %%

big_tab_list = []
sum_tot = defaultdict(int)
sum_active = defaultdict(int)
sum_act_agree = defaultdict(int)
sum_incon_agree = defaultdict(int)
sum_act_disagree = defaultdict(int)
sum_incon_disagree = defaultdict(int)

for clpath in cluster_paths:
    dset = clpath.parent.parent.stem

    df = pd.read_csv(clpath)

    act_tally = {}
    for _, r in pd.read_csv(switch_table_path).iterrows():
        if r["dataset"] == dset:
            act_tally[make_rowid(r)] = 0

    for _, r in df.iterrows():
        this = r["rowid"]

        mp = r["delta_mean_pval"]
        vp = r["delta_variance_pval"]

        del_var = r["delta_variance"]
        noise_del_var = r["noiseset_delta_variance"]

        if mp < mean_t and vp < var_t and del_var > noise_del_var:
            act_tally[this] += 1

    big_tab_list.extend(
        [
            {
                "microbe": dset,
                "rowid": k,
                "target_name": classname(k),
                "is_active": int(bool(v)),
                "lit_result": lit_agreement(k, int(bool(v))),
            }
            for k, v in act_tally.items()
        ]
    )

for d in big_tab_list:
    dset = d["microbe"]
    rid = d["rowid"]
    tname = d["target_name"]
    act = d["is_active"]
    lit = d["lit_result"]

    sum_tot[f"{dset}|{tname}"] += 1

    if lit is None:
        sum_active[f"{dset}|{tname}"] += act
        continue

    sum_act_agree[f"{dset}|{tname}"] += int(lit == "active_agree")
    sum_incon_agree[f"{dset}|{tname}"] += int(lit == "incon_agree")
    sum_act_disagree[f"{dset}|{tname}"] += int(lit == "active_disagree")
    sum_incon_disagree[f"{dset}|{tname}"] += int(lit == "incon_disagree")

big_tab = pd.DataFrame(big_tab_list)

# %%
# Merge literature expectations into the long inference table
big_tab = big_tab.merge(lit_table, how="outer", on=None)

# %%
big_tab.to_csv("./data/big_fig/locus_inferences.csv", index=False)

# %%
# Now make the summary table to build the circle plot
sum_list = []
for k, v in sum_tot.items():
    k: str
    s = k.split("|")
    sum_list.append(
        {
            "microbe": s[0],
            "target_name": s[1],
            "active": sum_active[k],
            "total": v,
            "sras": sra_dict[s[0]],
            "active_agree": sum_act_agree[k],
            "incon_agree": sum_incon_agree[k],
            "active_disagree": sum_act_disagree[k],
            "incon_disagree": sum_incon_disagree[k],
        }
    )

big_tab_sum = pd.DataFrame(sum_list)

# %%
big_tab_sum.to_csv("./data/big_fig/locus_inferences_sum.csv", index=False)

# %%
