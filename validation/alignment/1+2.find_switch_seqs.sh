#!/bin/bash

#SBATCH --time=00:15:00
#SBATCH --nodes=1 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/1+2.all.find_switch_seqs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 1+2.find_switch_seqs.sh

set echo on

RT="$HOME/packages/tart"
TABLE="$RT/validation/tables/inf_results.csv"
GENOMES="$RT/validation/genomes"

DSETS=(
    "a_baum"
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

PRE_DEL=1000
POST_DEL=1000

echo "Running reference generation and switch bounding for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${DSETS[@]}; do
    OUT_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta$PRE_DEL-$POST_DEL"

    mpiexec tart-targeted reference-gen \
        --ledger $TABLE \
        --out-dir $OUT_DIR \
        --genome $GENOMES \
        --dset $DSET \
        --unify \
        --pre-del $PRE_DEL \
        --post-del $POST_DEL &&
        wait

    echo "Finished reference generation for $DSET into $OUT_DIR"

    mpiexec tart-targeted bounds --ledger $TABLE \
        --out-dir $OUT_DIR \
        --genome $GENOMES \
        --dset $DSET \
        --pre-del $PRE_DEL \
        --post-del $POST_DEL &&
        wait

    echo "Finished switch bounds generation for $DSET into $OUT_DIR"
    echo
done

echo "Done"
