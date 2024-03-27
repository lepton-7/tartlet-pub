# Collects data for the big results figure from each peak plot cluster stats
#  file from the validation set

# %%
import pandas as pd

from glob import glob
from pathlib import Path


# %%
def classname(rowid: str):
    return rowid.split("#")[0]


mean_t = 0.01
var_t = 0.01

# %%
# Find all cluster stats files
cluster_paths = [
    Path(x) for x in glob("../alignment/outputs/*/plots/cluster_stats.csv")
]

# %%

big_tab_list = []
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

big_tab = pd.DataFrame(big_tab_list)
# %%
big_tab.to_csv("./data/big_fig/locus_inferences.csv", index=False)
# %%
