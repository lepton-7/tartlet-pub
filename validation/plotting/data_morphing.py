# %%
from typing import Any
import pandas as pd

from collections import defaultdict

# %%
# Dataset paths
cum_results_path = "../alignment/results.csv"

# %%
results_df = pd.read_csv(cum_results_path)

# %% --------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Make aliases for conditions and riboswitch loci
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
alias_dict = {}

for dset in pd.unique(results_df["dataset"]):
    subset = results_df[results_df["dataset"] == dset]
    prev = ""
    count = 1
    for rowid in pd.unique(subset["rowid"]):
        if prev != rowid.split("#")[0]:
            prev = rowid.split("#")[0]
            count = 1
        else:
            count += 1

        alias_dict[rowid] = f"{dset}_{prev}_{count}"

    prev = ""
    count = 1
    for cond in pd.unique(subset["condition"]):
        alias_dict[f"{dset}#{cond}"] = [dset, f"cond_{count}"]
        count += 1


# %% --------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Tally results by counting decisions
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

tally = defaultdict(lambda: defaultdict(dict))

for _, row in results_df.iterrows():
    id = row["rowid"]
    id_alias = alias_dict[row["rowid"]]
    dset = row["dataset"]
    cond = row["condition"]
    cond_alias = alias_dict[f"{dset}#{row['condition']}"][1]
    dec = row["decision"]
    geo = row["study"]

    try:
        tally[f"{id}"][f"{cond}"]["dataset"]
    except:
        tally[f"{id}"][f"{cond}"].update(
            {
                "geo": geo,
                "dataset": row["dataset"],
                "pass": 0,
                "fail": 0,
                "incon": 0,
            }
        )

    if dec != "pass" and dec != "fail":
        dec = "incon"

    # try:
    tally[f"{id}"][f"{cond}"][dec] += 1
    # except KeyError:

inter: list[dict[Any, Any]] = []

for rowid, conds in tally.items():
    for cond, decs in conds.items():
        t = {}
        t["switch_id"] = rowid
        t["condition"] = cond

        for decision, val in decs.items():
            t[decision] = val

        t["switch_alias"] = alias_dict[rowid]
        t["condition_alias"] = alias_dict[f"{decs['dataset']}#{cond}"][1]
        inter.append(t)

tally_df = pd.DataFrame.from_dict(inter)

# %%
tally_df.to_csv("./data/tally_aliased.csv", index=False)

# %%
