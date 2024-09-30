#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --nodes=3 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/5.all.make_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 5.make_BAMs.sh

set echo on

RT="$HOME/tartlet-pub"

DATECODE="20240930"

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

PRE_DEL=1000
POST_DEL=1000

echo "Making BAMs for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${DSETS[@]}; do

    REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta$PRE_DEL-$POST_DEL"

    SAM_DIR="$REF_DIR/alignment_$DATECODE"

    echo "Starting SAM conversion for $DSET"

    mpiexec tart-targeted convert-sam -i $SAM_DIR &&
        wait

    echo "Finished SAM -> sorted BAM conversion within $SAM_DIR"
    echo
done

echo "Done"
