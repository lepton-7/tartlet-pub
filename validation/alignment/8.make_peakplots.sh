#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=1 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/8.all.make_peakplots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 8.make_peakplots.sh

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

echo "Plotting output for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${DSETS[@]}; do

    D_ROOT="$RT/validation/alignment/outputs/$DSET"

    PPATH="$D_ROOT/plots/peak_log.csv"
    CPATH="$D_ROOT/plots/cluster_stats.csv"
    SAVE_PATH="$RT/validation/alignment/peak_plots/${DSET}_60.png"

    tartlet-targeted plot -p $PPATH \
        -c $CPATH \
        -o $SAVE_PATH \
        --name $DSET

    echo "----------------------------------------------------"
    echo
done

echo "Done"
