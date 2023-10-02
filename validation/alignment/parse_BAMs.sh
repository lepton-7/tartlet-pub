#!/bin/bash

#SBATCH --time=00:20:00
#SBATCH --nodes=5 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/parse_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j parse_BAMs.sh

set -x
set echo on

RT="$HOME/packages/tart"

DSET="b_theta"

BAM_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta500/alignments_2023-09-29_00-47-50"
SAVE_ROOT="$RT/validation/alignment/outputs/$DSET/plots/picks"
BOUNDS_DIR="$RT/validation/alignment/outputs/$DSET"

module load miniconda3
module load samtools
source activate local

mpiexec python parse_BAMs.py $BAM_DIR $SAVE_ROOT $BOUNDS_DIR &&
    wait
