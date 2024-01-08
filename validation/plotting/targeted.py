# %%
import pandas as pd

from matplotlib import pyplot as plt
from collections import defaultdict

# %%
# Dataset paths
cum_results_path = "../alignment/results.csv"

# %%
results_df = pd.read_csv(cum_results_path)

# %%
# Tally results by counting the decision column
tally = defaultdict(lambda: defaultdict(dict))

for _, row in results_df.iterrows():
    id = row["rowid"]
    cond = row["condition"]
    dec = row["decision"]
    tally[f"{id}"][f"{cond}"].update({"dataset": row["dataset"]})
    if dec != "pass" and dec != "fail":
        dec = "incon"

    try:
        tally[f"{id}"][f"{cond}"][dec] += 1
    except KeyError:
        tally[f"{id}"][f"{cond}"][dec] = 1

inter = []

for rowid, conds in tally.items():
    for cond, decs in conds.items():
        t = {}
        for decision, val in decs.items():
            t["switch_id"] = rowid
            t["condition"] = cond
            t[decision] = val
        inter.append(t)

tally_df = pd.DataFrame.from_dict(inter)

tally_df.to_csv("./data/tally.csv", index=False)