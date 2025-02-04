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
cl_sigpeak_dict = defaultdict(str)
both_pass_clsig_fail_dict = defaultdict(str)

for cpath in cluster_paths:
    df = pd.read_csv(cpath)

    for _, r in df.iterrows():
        rid = r["rowid"]

        mp = float(r["delta_mean_pval"])
        mval = float(r["delta_mean"])
        varp = float(r["delta_variance_pval"])
        varval = float(r["delta_variance"])
        noisevar = float(r["noiseset_delta_variance"])
        sigpeak_here = bool(r["sig_peak_in_cluster"])

        if mp < 0.05 and mval < 0:
            mean_dict[rid] = "pass"

        if varp < 0.05 and varval > noisevar:
            var_dict[rid] = "pass"

        if sigpeak_here:
            cl_sigpeak_dict[rid] = "pass"

        if (
            mp < 0.05
            and mval < 0
            and varp < 0.05
            and varval > noisevar
            and sigpeak_here
        ):
            both_dict[rid] = "pass"

        if (
            mp < 0.05
            and mval < 0
            and varp < 0.05
            and varval > noisevar
            and not sigpeak_here
        ):
            both_pass_clsig_fail_dict[rid] = "fail"

mean_dict_df = pd.DataFrame(mean_dict, index=["dec"]).T
var_dict_df = pd.DataFrame(var_dict, index=["dec"]).T
both_dict_df = pd.DataFrame(both_dict, index=["dec"]).T

sspeak_mean_passing_loci = [mid for mid in mean_dict.keys() if mid in passing_peak_loci]
sspeak_pass_mean_fail_loci = [
    rid for rid in passing_peak_loci if rid not in mean_dict.keys()
]
sspeak_mean_var_passing_loci = [
    vid for vid in var_dict.keys() if vid in sspeak_mean_passing_loci
]
sspeak_mean_pass_var_fail_loci = [
    rid for rid in sspeak_mean_passing_loci if rid not in var_dict.keys()
]
sspeak_mean_var_sig_passing_loci = [
    sigid for sigid in both_dict.keys() if sigid in sspeak_mean_var_passing_loci
]

print(f"Loci that have no passing peaks: {failcount}")
print(
    f"Loci that have passing peaks, but fail cl mean test: {len(sspeak_pass_mean_fail_loci)}"
)
print(
    f"Loci that have passing peaks, pass cl mean test, but fail cl var test: {len(sspeak_mean_pass_var_fail_loci)}"
)
print(
    f"Loci that fail despite passing cl mean+var test due to no passing peak in cluster: {len(both_pass_clsig_fail_dict)}"
)

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
# print(f"fail the mean fold coverage drop test: {395 - 18 - nacount - len(mean_dict)}")
# print(f"fail the variance test: {395 - 18 - nacount - len(var_dict)}")
mean_pass_var_fail = 0
mean_pass_var_pass = 0
passlist = []
for r in list(mean_dict.keys()):
    if var_dict[r] == "pass":
        mean_pass_var_pass += 1
        passlist.append(r)
    else:
        mean_pass_var_fail += 1
# print(f"passed mean test but failed variance test : {mean_pass_var_fail}")

# %%
# actiinf = infers[infers["is_active"] == 1]
# for i in passlist:
#     if i not in list(actiinf["rowid"]):
#         print(f"{i}")

# %%
# Check for which loci failed the locus plot decision but passed both cluster mean and variance tests
# for k in list(both_dict.keys()):
#     if k not in passing_peak_loci:
#         print(k)

# %%
sig_loc_list = [
    "Lysine#CP001139.1#1859200#1859017#-",
    "Lysine#CP080568.1#1148715#1148536#-",
    "Lysine#AE017334.2#1360986#1361167#+",
    "Lysine#AE017334.2#1696955#1696769#-",
    "Purine#AE017334.2#3605425#3605324#-",
    "SAM#AE017334.2#1375474#1375579#+",
    "SAM#AE017334.2#185584#185689#+",
    "SAM#AE017334.2#2953356#2953255#-",
    "SAM#AE017334.2#320605#320711#+",
    "SAM#AE017334.2#3890865#3890724#-",
    "SAM#AE017334.2#4740103#4739996#-",
    "SAM#AE017334.2#5140450#5140345#-",
    "TPP#AE017334.2#2433715#2433820#+",
    "TPP#AE017334.2#4953628#4953513#-",
    "TPP#AE017334.2#752099#752211#+",
    "glmS#AE017334.2#157527#157682#+",
    "ydaO-yuaA#AE017334.2#1768224#1768370#+",
    "ydaO-yuaA#AE017334.2#760272#760419#+",
    "ydaO-yuaA#AE017334.2#761408#761553#+",
    "yybP-ykoY#AE017334.2#1269290#1269126#-",
    "Cobalamin#CP069563.1#109168#108947#-",
    "Cobalamin#CP069563.1#1628761#1628968#+",
    "Cobalamin#CP069563.1#393861#393615#-",
    "Cobalamin#CP069563.1#4499045#4498827#-",
    "Cobalamin#CP008781.1#3811382#3811124#-",
    "Cobalamin#CP008781.1#733900#733634#-",
    "Glycine#CP008781.1#1833634#1833498#-",
    "Lysine#NC_000964.3#2911051#2910872#-",
    "Lysine#NC_000964.3#3421348#3421169#-",
    "M-box#NC_000964.3#1395616#1395781#+",
    "Purine#NC_000964.3#2320213#2320114#-",
    "Purine#NC_000964.3#698369#698470#+",
    "SAM#NC_000964.3#1180802#1180685#-",
    "SAM#NC_000964.3#1385891#1385736#-",
    "SAM#NC_000964.3#1424683#1424527#-",
    "SAM#NC_000964.3#1426876#1426976#+",
    "SAM#NC_000964.3#3129333#3129195#-",
    "SAM#NC_000964.3#3364503#3364397#-",
    "SAM#NC_000964.3#3997775#3997881#+",
    "SAM#NC_000964.3#3999166#3999272#+",
    "TPP#NC_000964.3#1607367#1607468#+",
    "ydaO-yuaA#NC_000964.3#486092#486235#+",
    "Cobalamin#CP040530.1#5410362#5410566#+",
    "Cobalamin#CP040530.1#5414425#5414630#+",
    "Cobalamin#CP040530.1#6034312#6034121#-",
    "Cobalamin#CP040530.1#6202554#6202342#-",
    "Cobalamin#CP072216.1#3081986#3081798#-",
    "Cobalamin#CP072216.1#3278312#3278102#-",
    "Cobalamin#CP076401.1#3260890#3260710#-",
    "Lysine#CP076401.1#1767336#1767512#+",
    "raiA#CP076401.1#185533#185740#+",
    "Cobalamin#CP001340.1#1950973#1951161#+",
    "Cobalamin#CP001340.1#528978#529205#+",
    "ZMP-ZTP#CP001340.1#88947#89027#+",
    "yybP-ykoY#CP001340.1#2597459#2597350#-",
    "Cobalamin#AE017285.1#721811#721631#-",
    "Cobalamin#NC_000913.3#4163384#4163573#+",
    "FMN#NC_000913.3#3184718#3184570#-",
    "MOCO_RNA_motif#NC_000913.3#816912#817056#+",
    "Mg_sensor#NC_000913.3#4467416#4467531#+",
    "TPP#NC_000913.3#2185432#2185336#-",
    "Lysine#KB944666.1#457062#456888#-",
    "M-box#KB944666.1#1484350#1484517#+",
    "Cobalamin#CP019962.1#4278947#4279124#+",
    "Cobalamin#CP019962.1#480456#480269#-",
    "MOCO_RNA_motif#CP019962.1#1479035#1479158#+",
    "Cobalamin#CP009494.1#6690206#6690014#-",
    "M-box#AL123456.3#2047595#2047769#+",
    "ydaO-yuaA#AL123456.3#965767#965545#-",
    "SAM#CP042324.1#4708468#4708615#+",
    "ydaO-yuaA#CP042324.1#3393037#3392871#-",
    "Cobalamin#AE006468.2#4347864#4348056#+",
    "FMN#AE006468.2#3358815#3358617#-",
    "FMN#CP035288.1#1086083#1086221#+",
    "Lysine#CP035288.1#1193293#1193468#+",
    "Lysine#CP035288.1#1493834#1493655#-",
    "AdoCbl_riboswitch#CP000387.1#512706#512558#-",
    "M-box#CP000387.1#1721468#1721312#-",
    "Cobalamin#CP046570.1#3179709#3179931#+",
]
