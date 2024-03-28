#!/bin/bash

#SBATCH --time=00:15:00
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
    "p_salmo"
    "p_cholor_aureo3084"
)

PRE_DEL=1000
POST_DEL=1000

echo "Running reference generation for:"
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
    echo
done

echo "Done"
