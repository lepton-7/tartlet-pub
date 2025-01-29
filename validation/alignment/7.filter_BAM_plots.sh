#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=5 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/7.all.filter_BAM_plots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 7.filter_BAM_plots.sh

set echo on

RT="$HOME/tartlet-pub"

DSETS=(
    "a_fischeri_ES114"
    "a_kunk"
    "b_anth"
    "b_frag"
    "b_pseudo"
    "b_sub_168"
    "b_theta"
    "b_xyla"
    "c_diff"
    "c_vibrioides"
    "d_vulg"
    "e_coli"
    "e_fae"
    "e_limo"
    "k_pneum"
    "m_smeg"
    "m_tuber"
    "n_gonorr"
    "p_cholor_aureo3084"
    "p_fluor"
    "p_salmo"
    "s_coelicolor"
    "s_elon"
    "s_enter_typh"
    "s_epi"
    "s_meli"
    "s_sanguinis"
    "s_spcc6803"
    "x_albi"
    "x_ory"
)
# These were removed from ^
#   "c_basil"
#   "p_aeru"

echo "Filtering output plots for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${DSETS[@]}; do

    D_ROOT="$RT/validation/alignment/outputs/$DSET"

    PICKLE_ROOT="$D_ROOT/plots/picks.tar.gz"
    SAVE_ROOT="$D_ROOT/plots/"

    # echo "Removing existing data"
    # rm -r $SAVE_ROOT/*/ $SAVE_ROOT/*.csv
    rm -r $SAVE_ROOT/*.csv

    echo "Filtering plots for $DSET"

    mpiexec tartlet-targeted filter -i $PICKLE_ROOT \
        -o $SAVE_ROOT \
        --ext-prop -0.3 1.0 \
        --conv \
        --min-cov-depth 15 \
        --cophen-dist-thresh 0.04 \
        --noplots \
    # --statplot \

    echo "Finished filtering outputs for $DSET into $SAVE_ROOT"
    echo "----------------------------------------------------"
    echo
done

echo "Done"
