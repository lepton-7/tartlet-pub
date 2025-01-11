# %%
import pandas as pd

from pathlib import Path
from collections import defaultdict

rt = "../.."

# %%
ppaths = [*Path(f"{rt}/validation/alignment/outputs").glob("**/plots/peak_log.csv")]

sig_clust_map = defaultdict(int)
tot_pass_peaks = 0
tot_pass_in_sig_clust = 0

# %%
for ppath in ppaths:
    cpath = ppath.parent.joinpath("cluster_stats.csv")

    p = pd.read_csv(ppath)
    c = pd.read_csv(cpath)

    c_sub = c[
        (c["delta_mean_pval"] < 0.05)
        & (c["delta_variance_pval"] < 0.05)
        & (c["delta_variance"] > c["noiseset_delta_variance"])
    ]

    for _, r in c_sub.iterrows():
        sig_clust_map[f"{r['rowid']}#{r['cluster']}"] = 1

    for _, r in p[p["decision"] == "pass"].iterrows():
        tot_pass_peaks += 1
        tot_pass_in_sig_clust += sig_clust_map[f"{r['rowid']}#{r['cluster']}"]

print(f"{tot_pass_peaks} passing peaks; {tot_pass_in_sig_clust} in sig clusters")
