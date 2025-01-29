# %%
from typing import Literal
import pandas as pd

from glob import glob
from pathlib import Path
from collections import defaultdict


# %%
excluded_dsets = [
    "a_baum",
    "c_basil",
    "p_fluor",
    "p_aeru",
    "p_salmonis",
]

peak_paths = [
    Path(x)
    for x in glob("./outputs/*/plots/peak_log.csv")
    if Path(x).parent.parent.stem not in excluded_dsets
]
cluster_paths = [
    Path(x)
    for x in glob("./outputs/*/plots/cluster_stats.csv")
    if Path(x).parent.parent.stem not in excluded_dsets
]

# %%

peak_dict = defaultdict(str)
for ppath in peak_paths:
    df = pd.read_csv(ppath)

    for _, r in df.iterrows():
        rid = r["rowid"]

        if peak_dict[rid] == "pass":
            continue
        peak_dict[rid] = r["decision"]

peak_dec_df = pd.DataFrame(peak_dict, index=["dec"]).T
passing_peak_loci = list(peak_dec_df[peak_dec_df["dec"] == "pass"].index)

pcount = 0  # Pass decision counter
nacount = 0  # This indicates how many failed because of insufficient coverage
failcount = 0  # Fail decision counter

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
both_dict = defaultdict(str)
for cpath in cluster_paths:
    df = pd.read_csv(cpath)

    for _, r in df.iterrows():
        rid = r["rowid"]

        mp = float(r["delta_mean_pval"])
        mval = float(r["delta_mean"])
        varp = float(r["delta_variance_pval"])
        varval = float(r["delta_variance"])
        noisevar = float(r["noiseset_delta_variance"])

        if mp < 0.05 and mval < 0:
            mean_dict[rid] = "pass"

        if varp < 0.05 and varval > noisevar:
            var_dict[rid] = "pass"

        if mp < 0.05 and mval < 0 and varp < 0.05 and varval > noisevar:
            both_dict[rid] = "pass"


# %%
# Test why there are fewer loci in the peak_dict than the expected 395

rt = "../.."
infers = pd.read_csv(f"{rt}/validation/plotting/data/big_fig/locus_inferences.csv")
undercount = 0
for _, r in infers.iterrows():
    k = r["rowid"]
    if k not in list(peak_dict.keys()):
        print(f"{r['microbe']}: {k}")
        undercount += 1

# %%
print(f"fail the mean fold coverage drop test: {395 - 18 - nacount - len(mean_dict)}")
print(f"fail the variance test: {395 - 18 - nacount - len(var_dict)}")
mean_pass_var_fail = 0
mean_pass_var_pass = 0
passlist = []
for r in list(mean_dict.keys()):
    if var_dict[r] == "pass":
        mean_pass_var_pass += 1
        passlist.append(r)
    else:
        mean_pass_var_fail += 1
print(f"passed mean test but failed variance test : {mean_pass_var_fail}")

# %%
actiinf = infers[infers["is_active"] == 1]
for i in passlist:
    if i not in list(actiinf["rowid"]):
        print(f"{i}")

# %%
# Check for which loci failed the locus plot decision but passed both cluster mean and variance tests
for k in list(both_dict.keys()):
    if k not in passing_peak_loci:
        print(k)
