# %%
import pandas as pd

from pathlib import Path
from collections import defaultdict

rt = "../.."

# %%
ppaths_whole = [
    *Path(f"{rt}/validation/alignment/outputs").glob("**/plots/peak_log.csv")
]

excluded_dsets = [
    "a_baum",
    "c_basil",
    "p_fluor",
    "p_aeru",
    "p_salmonis",
]

ppaths = [x for x in ppaths_whole if x.parent.parent.stem not in excluded_dsets]


# %%

sig_clust_map = defaultdict(int)
tot_pass_peaks = 0
tot_pass_in_sig_clust = 0
tot_peaks = 0


for ppath in ppaths:
    cpath = ppath.parent.joinpath("cluster_stats.csv")

    p = pd.read_csv(ppath)
    c = pd.read_csv(cpath)

    tot_peaks += p.shape[0]

    c_sub = c[
        (c["delta_mean_pval"] < 0.05)
        & (c["delta_variance_pval"] < 0.05)
        & (c["delta_mean"] < 0)
        & (c["delta_variance"] > c["noiseset_delta_variance"])
        & (c["sig_peak_in_cluster"] == True)
    ]

    for _, r in c_sub.iterrows():
        sig_clust_map[f"{r['rowid']}#{r['cluster']}"] = 1

    for _, r in p[p["decision"] == "pass"].iterrows():
        tot_pass_peaks += 1
        tot_pass_in_sig_clust += sig_clust_map[f"{r['rowid']}#{r['cluster']}"]

print(
    f"{tot_peaks} peaks in total; {tot_pass_peaks} passing peaks; {tot_pass_in_sig_clust} in {sum(sig_clust_map.values())} sig clusters"
)

# %%

sigrid_list = [
    "#".join(k.split("#")[:-1]) for k in sig_clust_map.keys() if sig_clust_map[k] > 0
]
sigrid_list = list(pd.unique(sigrid_list))

# %%
