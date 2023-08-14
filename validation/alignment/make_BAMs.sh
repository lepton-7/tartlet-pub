#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=10 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/make_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j make_BAMs.sh

set -x
set echo on

RT="$HOME/packages/tart"
REF_DIR="$RT/validation/alignment/switch_seqs_delta500"
SAM_DIR="$REF_DIR/alignments_2023-08-14_11-02-43"

module load miniconda3
module load samtools
source activate local

mpiexec python -u make_BAMs.py $SAM_DIR &&
    wait
