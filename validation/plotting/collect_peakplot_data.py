# Collects data for the big results figure from each peak plot cluster stats
#  file from the validation set

# %%
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict


# %%
def classname(rowid: str):
    return rowid.split("#")[0]


mean_t = 0.01
var_t = 0.01

# %%
# Find all cluster stats files
cluster_paths = [
    # Path(x) for x in glob("../alignment/outputs/*/plots/cluster_stats.csv")
    Path(x) for x in glob("../alignment/outputs/*/p2/cluster_stats.csv")
]

# %%

big_tab_list = []
sum_tot = defaultdict(int)
sum_active = defaultdict(int)
for clpath in cluster_paths:
    dset = clpath.parent.parent.stem

    df = pd.read_csv(clpath)

    act_tally = {}
    for rid in pd.unique(df["rowid"]):
        act_tally[rid] = 0

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
            }
            for k, v in act_tally.items()
        ]
    )
    for k, v in act_tally.items():
        sum_tot[f"{dset}|{classname(k)}"] += 1
        sum_active[f"{dset}|{classname(k)}"] += int(bool(v))

big_tab = pd.DataFrame(big_tab_list)
# %%
big_tab.to_csv("./data/big_fig/locus_inferences.csv", index=False)

# %%
# Now make the summary table to build the circle plot
sum_list = []
for k, v in sum_tot.items():
    k: str
    s = k.split("|")
    sum_list.append(
        {"microbe": s[0], "target_name": s[1], "active": sum_active[k], "total": v}
    )

big_tab_sum = pd.DataFrame(sum_list)

# %%
big_tab_sum.to_csv("./data/big_fig/locus_inferences_sum.csv", index=False)
