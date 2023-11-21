#!/bin/bash

#SBATCH --time=00:50:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/2.all.switch_loc_in_ref.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 2.switch_loc_in_ref.sh

set echo on

RT="$HOME/packages/tart"
TABLE="$RT/validation/tables/inf_results.csv"
GENOMES="$RT/validation/genomes"

DSETS=(
    "c_basil"
    "m_smeg"
    "s_meli"
    "d_vulg"
    "b_frag"
    "e_fae"
)

PRE_DEL=1000
POST_DEL=1000

echo "Finding switch locations in reference for:"
echo "${DSETS[*]}"
echo
echo

for DSET in ${DSETS[@]}; do
    OUT_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta$PRE_DEL-$POST_DEL"

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
