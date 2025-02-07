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
pdf = pd.DataFrame()
for ppath in peak_paths:
    df = pd.read_csv(ppath)
    pdf = pd.concat([pdf, df])
    for _, r in df.iterrows():
        rid = r["rowid"]

        if peak_dict[rid] == "pass":
            continue

        peak_dict[rid] = r["decision"]

pdf_nonaloc = pdf[pdf["decision"].isin(["fail", "pass"])]
nancount = len(pd.unique(pdf["rowid"])) - len(pd.unique(pdf_nonaloc["rowid"]))

pdf_passingloc = pdf_nonaloc[pdf_nonaloc["decision"] == "pass"]
pcount = len(pd.unique(pdf_passingloc["rowid"]))

failcount = len(pd.unique(pdf_nonaloc["rowid"])) - pcount

peak_dec_df = pd.DataFrame(peak_dict, index=["dec"]).T
passing_peak_loci = list(peak_dec_df[peak_dec_df["dec"] == "pass"].index)
failing_peak_loci = list(peak_dec_df[peak_dec_df["dec"] == "fail"].index)

# %%

mean_dict = defaultdict(str)
var_dict = defaultdict(str)
both_dict = defaultdict(str)
cl_sigpeak_dict = defaultdict(str)
both_pass_clsig_fail_dict = defaultdict(str)
cldf = pd.DataFrame()

for cpath in cluster_paths:
    df = pd.read_csv(cpath)
    cldf = pd.concat([cldf, df])

cldf["sigcl"] = (
    (cldf["delta_mean_pval"] < 0.05)
    & (cldf["delta_mean"] < 0)
    & (cldf["delta_variance"] > cldf["noiseset_delta_variance"])
    & (cldf["sig_peak_in_cluster"])
    & (cldf["delta_variance_pval"] < 0.05)
)

cldf["meanpass"] = (cldf["delta_mean_pval"] < 0.05) & (cldf["delta_mean"] < 0)

cldf["varpass"] = (cldf["delta_variance"] > cldf["noiseset_delta_variance"]) & (
    cldf["delta_variance_pval"] < 0.05
)

tot_rowid_in_cstats = len(pd.unique(cldf["rowid"]))
print(f"Total loci represented across all cstats: {tot_rowid_in_cstats}")


ppas_rowids = list(pd.unique(cldf[cldf["sig_peak_in_cluster"]]["rowid"]))
cldf_ppassing_loci = cldf[cldf["rowid"].isin(ppas_rowids)]
ppass_rowids_count = len(pd.unique(cldf_ppassing_loci["rowid"]))
print(f"Loci that have no sig peaks: {tot_rowid_in_cstats-ppass_rowids_count}\n")
print(
    f"Therefore, loci that have at least one sig peak: {ppass_rowids_count}\nGating on this condition:\n"
)


meandf = cldf_ppassing_loci[
    cldf_ppassing_loci["meanpass"] & cldf_ppassing_loci["sig_peak_in_cluster"]
]
meanpassing_loc = len(pd.unique(meandf["rowid"]))
print(f"Loci with atleast one sig-peak cl passing mean test: {meanpassing_loc}")
print(
    f"Therefore, loci that have no sig-peak cls that are mean passing: {ppass_rowids_count - meanpassing_loc}\n"
)

mean_var_df = meandf[meandf["varpass"]]
meanvarpassing_loc = len(pd.unique(mean_var_df["rowid"]))
print(
    f"Loci with atleast one sig-peak cl passing mean and var test: {meanvarpassing_loc}"
)
print(
    f"Therefore, loci that have no sig-peak cls that pass both mean+var: {meanpassing_loc - meanvarpassing_loc}\n"
)

# %%
