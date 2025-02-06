#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=1 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/8.all.make_peakplots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 8.make_peakplots.sh

set echo on

RT="$HOME/tartlet-pub"

declare -A DSETS=(
    ["a_fischeri_ES114"]="A. fischeri"
    ["a_kunk"]="A. kunkeei"
    ["b_anth"]="B. anthracis"
    ["b_frag"]="B. fragilis"
    ["b_pseudo"]="B. pseudomallei"
    ["b_sub_168"]="B. subtilis"
    ["b_theta"]="B. thetaiotaomicron"
    ["b_xyla"]="B. xylanisolvens"
    ["c_diff"]="C. difficile"
    ["c_vibrioides"]="C. vibrioides"
    ["d_vulg"]="D. vulgaris"
    ["e_coli"]="E. coli"
    ["e_fae"]="E. faecalis"
    ["e_limo"]="E. limosum"
    ["k_pneum"]="K. pneumoniae"
    ["m_smeg"]="M. smegmatis"
    ["m_tuber"]="M. tuberculosis"
    ["n_gonorr"]="N. gonorrhoeae"
    ["p_cholor_aureo3084"]="P. chlororaphis"
    ["p_fluor"]="P. fluorescens"
    ["p_salmo"]="P. salmonis"
    ["s_coelicolor"]="S. coelicolor"
    ["s_elon"]="S. elongatus"
    ["s_enter_typh"]="S. enterica"
    ["s_epi"]="S. epidermidis"
    ["s_meli"]="S. meliloti"
    ["s_sanguinis"]="S. sanguinis"
    ["s_spcc6803"]="Synechococcus sp."
    ["x_albi"]="X. albilineans"
    ["x_ory"]="X. oryzae"
)
# These were removed from ^
#   "c_basil"
#   "p_aeru"

echo "Plotting output for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${!DSETS[@]}; do

    D_ROOT="$RT/validation/alignment/outputs/$DSET"

    PPATH="$D_ROOT/plots/peak_log.csv"
    CPATH="$D_ROOT/plots/cluster_stats.csv"
    SAVE_PATH="$RT/validation/alignment/peak_plots/${DSET}.png"

    tartlet-targeted plot -p $PPATH \
        -c $CPATH \
        -o $SAVE_PATH \
        --name "${DSETS[$DSET]}"

    echo "----------------------------------------------------"
    echo
done

echo "Done"
