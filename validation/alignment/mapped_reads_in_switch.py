# %%
from pathlib import Path
from io import TextIOWrapper
from collections import defaultdict

# %%
rt = "../.."
file_list = [
    Path(
        f"{rt}/validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-48-10.31786767"
    ),
    Path(
        f"{rt}/validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-17.31786806"
    ),
    Path(
        f"{rt}/validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-29.31786821"
    ),
    Path(
        f"{rt}/validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-43.31786822"
    ),
]


def get_mapped_rates(f: TextIOWrapper):
    d = defaultdict(list)
    cur_set = ""
    for l in f:
        if l.startswith("Starting BAM parsing"):
            cur_set = l.split(" ")[-1].rstrip()
            d[cur_set] = [0, 0]
            print(f"Calculating mapped rates for {cur_set}")

        if l.startswith("Mapped in switch"):
            mline = l.split(";")[0]
            uline = l.split(";")[1]
            m = int(mline.split(" = ")[-1])
            u = int(uline.split(" = ")[-1])

            d[cur_set][0] += m
            d[cur_set][1] += u

        if l.startswith("Finished BAM parsing"):
            print(
                f"{d[cur_set][0]} mapped and {d[cur_set][1]} unmapped in switch across {cur_set}\n"
            )

    return d


# %%
rates_dict = {}
for filepath in file_list:
    with open(filepath, "r") as f:
        rates_dict.update(get_mapped_rates(f))

# %%
total_mapped = 0
total_unmapped = 0

for k, vlist in rates_dict.items():
    total_mapped += vlist[0]
    total_unmapped += vlist[1]

map_rate = total_mapped / (total_mapped + total_unmapped)

print(
    f"In total: {total_mapped:,} mapped and {total_unmapped:,} unmapped ({(map_rate*100):0.2f}%)"
)
