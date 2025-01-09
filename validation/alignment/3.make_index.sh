#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/3.all.make_index.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 3.make_index.sh

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
    "c_basil"
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
    "p_aeru"
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

for DSET in ${DSETS[@]}; do

    REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000"

    srun -n 1 tartlet-targeted index -i $REF_DIR -p 44

    echo "Finished index generation for $DSET from $REF_DIR"
    echo
done
