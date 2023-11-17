#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=3 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/5.m_tuber.make_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 5.make_BAMs.sh

set echo on

RT="$HOME/packages/tart"

DSET="m_tuber"

REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000"

SAM_DIR="$REF_DIR/alignment_20231113"

echo "Starting SAM conversion for $DSET"

mpiexec tart-targeted convert-sam -i $SAM_DIR &&
    wait

echo "Finished SAM -> sorted BAM conversion within $SAM_DIR"
