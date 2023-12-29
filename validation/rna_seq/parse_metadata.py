# %%
import pandas as pd

from glob import glob
from pathlib import Path
from shutil import rmtree

# %%
keep_cols = {
    "GSE101911": ["genotype"],
    "GSE114358": ["carbon_source"],
    "GSE122211": ["knock_out", "media"],
    "GSE122295": ["supplement", "media"],
    "GSE122296": ["strain", "media"],
    "GSE122779": ["treatment", "media"],
    "GSE129572": ["carbon_source"],
    "GSE132275": ["source_name", "light_period"],
    "GSE149269": ["Substrate", "growth_phase"],
    "GSE151705": ["treatment", "strain"],
    "GSE152295": ["treatment"],
    "GSE152356": ["genotype"],
    "GSE169260": ["genotype", "dp_treatment"],
    "GSE169409": ["growth_phase", "strain"],
    "GSE173804": ["atmosphere"],
    "GSE174672": ["strain", "Timepoint"],
    "GSE183334": ["genotype"],
    "GSE191020": ["Substrate", "treatment"],
    "GSE200822": ["treatment", "Time"],
    "GSE203342": ["treatment", "genotype"],
    "GSE205998": ["source_name"],
    "GSE218354": ["treatment", "vector"],
    "GSE219221": ["treatment"],
    "GSE220692": ["treatment"],
    "GSE227397": ["treatment", "genotype"],
    "GSE229478": ["genotype"],
    "GSE229867": ["strain"],
    "GSE232901": ["stress", "Time_point"],
    "GSE234439": ["treatment", "genotype"],
    "GSE235725": ["strain"],
    "GSE237189": ["genotype", "treatment"],
    "GSE241057": ["genotype", "treatment"],
    "GSE61200": ["strain"],
    "GSE74379": ["carbon_source", "growth_phase"],
    "GSE78834": ["fna_level"],
    "GSE89651": ["source_name"],
}
# %%
m_list = glob("metadata/*.txt")

for x in m_list:
    df = pd.read_csv(x)
    df.to_csv(f"{x[:-4]}.csv", index=False)
    # rmtree(Path(x))
# %%
conditions = {}
for csv in glob("metadata/*.csv"):
    name = Path(csv).stem
    df = pd.read_csv(csv)

    cols = keep_cols[name]
    for _, row in df.iterrows():
        conditions[row["Run"]] = (";".join([str(row[x]) for x in cols]), name)

# %%
df = pd.DataFrame.from_dict(conditions, orient="index", columns=["condition", "study"])
df.index.names = ["Run"]
df.to_csv("run_conditions.csv")
