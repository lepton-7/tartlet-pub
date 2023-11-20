#!/bin/bash

#SBATCH --time=00:50:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/1.all.find_switch_seqs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 1.find_switch_seqs.sh

set echo on

RT="$HOME/packages/tart"
TABLE="$RT/validation/tables/inf_results.csv"
GENOMES="$RT/validation/genomes"

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

echo "Running reference generation for:"
echo "${DSETS[*]}"

for DSET in ${DSETS[@]}; do
    DSET="a_rabiei"

    PRE_DEL=1000
    POST_DEL=1000

    OUT_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta$PRE_DEL-$POST_DEL"

    mpiexec tart-targeted reference-gen --ledger $TABLE \
        --out-dir $OUT_DIR \
        --genome $GENOMES \
        --dset $DSET \
        --pre-del $PRE_DEL \
        --post-del $POST_DEL &&
        wait

    echo "Finished reference generation for $DSET into $OUT_DIR"
done
