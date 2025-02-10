# %%
from pathlib import Path
from io import TextIOWrapper
from collections import defaultdict

# %%
rt = "../.."
unsoft_list = [
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-48-10.31786767",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-17.31786806",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-29.31786821",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-09-30_10-50-43.31786822",
]
soft_list = [
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-10-01_00-37-04.31801983",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-10-01_00-37-20.31801984",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-10-01_00-37-33.31801985",
    "validation/alignment/jobs/6.all.parse_BAMs.out.2024-10-01_00-37-55.31801991",
]

file_list = [Path(f"{rt}/{relpath}") for relpath in soft_list]


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
excluded_dsets = [
    "a_baum",
    "c_basil",
    "p_fluor",
    "p_aeru",
    "p_salmo",
]

for k, vlist in rates_dict.items():
    if k in excluded_dsets:
        print(f"Excluding {k} from total counts")
        continue
    total_mapped += vlist[0]
    total_unmapped += vlist[1]

map_rate = total_mapped / (total_mapped + total_unmapped)

print(
    f"In total: {total_mapped:,} mapped and {total_unmapped:,} unmapped ({(map_rate*100):0.2f}%)"
)
