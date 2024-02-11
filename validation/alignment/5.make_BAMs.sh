#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=3 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/5.all.make_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 5.make_BAMs.sh

set echo on

RT="$HOME/packages/tart"

# DATECODE="20231113"
DATECODE="20240209"

DSETS=(
    "b_sub_168"
    # "n_gonorr"
)

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
