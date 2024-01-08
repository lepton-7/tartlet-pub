# %%
import pandas as pd

from pandas.errors import EmptyDataError
from glob import glob
from pathlib import Path

# %%
res_paths = glob("outputs/*/plots/characteristics.*")
cond_path = "../rna_seq/run_conditions.csv"

cond_df = pd.read_csv(cond_path)
# %%
# Compile output characteristics into one file

df = pd.DataFrame()
for dfpath in res_paths:
    try:
        d2 = pd.read_csv(dfpath)
    except EmptyDataError:
        print(f"No data in {dfpath}")
        continue
    bug = Path(dfpath).parent.parent.stem
    buglist = [bug for i in range(len(d2))]
    d2["dataset"] = buglist

    df = pd.concat([df, d2], ignore_index=True)

# %%
merged = df.merge(
    cond_df,
    left_on="transcriptome",
    right_on="Run",
    right_index=False,
    left_index=False,
)
merged.drop(columns=["Run"], inplace=True)

# %%
merged.to_csv("results.csv", index=False)
