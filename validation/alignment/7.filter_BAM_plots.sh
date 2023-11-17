#!/bin/bash

#SBATCH --time=00:50:00
#SBATCH --nodes=3 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/7.all.filter_BAM_plots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 7.filter_BAM_plots.sh

# set -x
set echo on

RT="$HOME/packages/tart"

DSETS=(
    "a_baum"
    "a_rabiei"
    "c_vibrioides"
    "m_tuber"
    "s_elon"
    "s_spcc6803"
    "a_fischeri_ES114"
    "b_sub_168"
    "e_coli"
    "p_fluor"
    "s_enter_typh"
    "x_albi"
    "a_kunk"
    "b_theta"
    "k_pneum"
    "s_coelicolor"
    "s_sanguinis"
)

for DSET in ${DSETS[@]}; do

    D_ROOT="$RT/validation/alignment/outputs/$DSET"

    PICKLE_ROOT="$D_ROOT/plots/picks.tar.gz"
    SAVE_ROOT="$D_ROOT/plots/"

    echo "Filtering plots for $DSET"

    mpiexec tart-targeted filter -i $PICKLE_ROOT -o $SAVE_ROOT --ext-prop -0.3 1.0 --conv --min-cov-depth 25

    echo "Finished filtering outputs for $DSET into $SAVE_ROOT"
    echo "----------------------------------------------------"
    echo
done
