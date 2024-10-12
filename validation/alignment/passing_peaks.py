# %%
from typing import Literal
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict


# %%
peak_paths = [Path(x) for x in glob("./outputs/*/plots/peak_log.csv")]
cluster_paths = [Path(x) for x in glob("./outputs/*/plots/cluster_stats.csv")]
# %%

peak_dict = defaultdict(str)
for ppath in peak_paths:
    df = pd.read_csv(ppath)

    for _, r in df.iterrows():
        rid = r["rowid"]

        if peak_dict[rid] == "pass":
            continue
        peak_dict[rid] = r["decision"]

pcount = 0
nacount = 0
failcount = 0

for k, v in peak_dict.items():
    if v == "fail":
        failcount += 1
    elif v == "pass":
        pcount += 1
    else:
        nacount += 1

# %%

mean_dict = defaultdict(str)
var_dict = defaultdict(str)
for cpath in cluster_paths:
    df = pd.read_csv(cpath)

    for _, r in df.iterrows():
        rid = r["rowid"]

        mp = float(r["delta_mean_pval"])
        mval = float(r["delta_mean"])
        varp = float(r["delta_variance_pval"])

        if mp < 0.05 and mval < 0:
            mean_dict[rid] = "pass"

        if varp < 0.05:
            var_dict[rid] = "pass"
