# %%

from pathlib import Path
from tartlet.targeted.plot_results import main
from tartlet.utils.mpi_context import BasicMPIContext

rt = "../.."
# %%
cophens = [x / 100 for x in range(1, 30)]

mp_con = BasicMPIContext()
comm = mp_con.comm
rank = mp_con.rank

dset = "b_sub_168"

# %%
for cophen in cophens:
    ppath = Path(f"outputs_{dset}/cophen{cophen:.2f}/peak_log.csv")
    cpath = Path(f"outputs_{dset}/cophen{cophen:.2f}/cluster_stats.csv")

    outpath = Path(f"plots/{dset}_{cophen:.2f}.png")

    main(str(ppath), str(cpath), f"{dset}_d={cophen:.2f}", str(outpath))

    if rank == 0:
        print(f"Generated {dset} peak plots for clustering distance threshold = {cophen:.2f}")

# %%
