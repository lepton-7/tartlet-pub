# %%

from pathlib import Path
from tartlet.targeted.filter_BAM_plots import main
from tartlet.utils.mpi_context import BasicMPIContext

rt = "../.."
# %%
cophens = [x / 100 for x in range(1, 10)]

mp_con = BasicMPIContext()
comm = mp_con.comm
rank = mp_con.rank

dset = "e_coli"

# %%
for cophen in cophens:
    odir = Path(f"outputs_{dset}/cophen{cophen:.2f}")
    odir.mkdir(parents=True, exist_ok=True)

    main(
        pick_root=f"{rt}/validation/alignment/outputs/{dset}/plots/picks.tar.gz",
        out_dir=odir,
        ext_prop=(-0.3, 1.0),
        roi_val=0.5,
        noplots=True,
        cophen_dist_thresh=cophen,
        bin_size=10,
        min_cov_depth=15,
        conv=True,
        statplot=False,
        rel_cov_change_sig_thresh=-0.2,
    )
    if rank == 0:
        print(
            f"Generated outputs for {dset} clustering distance threshold = {cophen:.2f}"
        )

# %%
