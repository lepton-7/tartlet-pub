# %%
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.stats import norm
from matplotlib.axes import Axes
from matplotlib.text import Text

# %%
palette = {
    "Read": "#3D5A80",
    "axback": "#FFFFFF",
    "figback": "#FFFFFF",
}

# %%
# Data for plotting

arr_size = 70
x_arr = [i for i in range(arr_size)]

raw = np.zeros(arr_size)
raw[16] = 1
raw[17] = -1
raw[60] = 1

num = 51 // 2
k_in = [x for x in range(-num, num + 1)]
kernel = norm.pdf(k_in, loc=0, scale=3)

conv = np.convolve(raw, kernel, "same")

# %%
fontsize = "xx-large"


def gen_formatting(ax: Axes):
    ax.set_xlabel("Position along reference", fontsize=fontsize)

    xticklabs = ["p", "p+1", "t"]
    ax.set_xticks([16, 17, 60])
    ax.set_xticklabels(xticklabs, fontsize="large")
    tlabs: list[Text] = ax.get_xticklabels()
    tlabs[0].set_ha("right")  # type: ignore
    tlabs[1].set_ha("left")  # type: ignore

    ax.set_yticks([])

    ax.axhline(linewidth=0.8, y=0, color="black")

    ax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# %%

fig, ax = plt.subplots(
    1,
    2,
    sharex=True,
    figsize=(10, 3),
    dpi=330,
    constrained_layout=True,
    facecolor=palette["figback"],
)

rawplt: Axes = ax[0]
convplt: Axes = ax[1]

rawplt.bar(
    x_arr,
    raw,
    width=1,
    color=palette["Read"],
    align="edge",
)
rawplt.set_ylabel("Sum of termini", fontsize=fontsize)
gen_formatting(rawplt)

convplt.bar(
    x_arr,
    conv,
    width=1,
    color=palette["Read"],
    align="edge",
)
convplt.set_ylabel("Convolved signed \nsum", fontsize=fontsize)
convplt.set_ylim((-0.14, 0.14))
gen_formatting(convplt)


# %%

save_path = "./plots/concept_figs/__conv_application.svg"
fig.savefig(save_path)
