#!/bin/bash

#SBATCH --time=00:15:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/filter_BAM_plots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j filter_BAM_plots.sh

set -x
set echo on

RT="$HOME/packages/tart"

DSET="p_fluor"

D_ROOT="$RT/validation/alignment/outputs/$DSET"

PICKLE_ROOT="$D_ROOT/plots/picks"
SAVE_ROOT="$D_ROOT/plots/filtered"

module load miniconda3
source activate local

mpiexec python filter_BAM_plots.py $PICKLE_ROOT $SAVE_ROOT &&
    wait
