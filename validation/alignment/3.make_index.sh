#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/make_index.out.$(date +"%Y-%m-%d_%H-%M-%S").%j make_index.sh

set echo on

DSETS=("a_fischeri_ES114" "a_rabiei" "b_theta" "c_vibrioides" "s_coelicolor" "s_spcc6803")

for DSET in ${DSETS[@]}; do

    RT="$HOME/packages/tart"
    REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta500"

    srun -n 1 tart-targeted index -i $REF_DIR -p 44

    echo "Finished index generation for $DSET from $REF_DIR"
    echo
done
