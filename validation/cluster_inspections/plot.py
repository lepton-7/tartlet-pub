# %%

from pathlib import Path
from tartlet.targeted.plot_results import main
from tartlet.utils.mpi_context import BasicMPIContext

rt = "../.."
# %%
cophens = [x / 100 for x in range(1, 100)]

mp_con = BasicMPIContext()
comm = mp_con.comm
rank = mp_con.rank

# %%
for cophen in cophens:
    ppath = Path(f"outputs/cophen{cophen}/peak_log.csv")
    cpath = Path(f"outputs/cophen{cophen}/cluster_stats.csv")

    outpath = Path(f"plots/e_coli_{cophen}.png")

    main(str(ppath), str(cpath), f"e_coli_d={cophen}", str(outpath))

    if rank == 0:
        print(f"Generated peak plots for clustering distance threshold = {cophen}")

# %%
