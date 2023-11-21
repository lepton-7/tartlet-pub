#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/3.all.make_index.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 3.make_index.sh

set echo on

DSETS=(
    "c_basil"
    "m_smeg"
    "s_meli"
    "d_vulg"
    "b_frag"
    "e_fae"
)


for DSET in ${DSETS[@]}; do

    RT="$HOME/packages/tart"
    REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000"

    srun -n 1 tart-targeted index -i $REF_DIR -p 44

    echo "Finished index generation for $DSET from $REF_DIR"
    echo
done
