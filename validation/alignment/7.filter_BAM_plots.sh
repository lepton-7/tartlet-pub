#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/filter_BAM_plots.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 7.filter_BAM_plots.sh

# set -x
set echo on

RT="$HOME/packages/tart"
DSET="b_sub_168"

D_ROOT="$RT/validation/alignment/outputs/$DSET"

PICKLE_ROOT="$D_ROOT/plots/picks"
SAVE_ROOT="$D_ROOT/plots/"

mpiexec tart-targeted filter -i $PICKLE_ROOT -o $SAVE_ROOT --ext-prop -0.3 1.0 --conv

echo "Finished filtering outputs for $DSET into $SAVE_ROOT"
