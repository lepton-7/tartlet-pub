#!/bin/bash

#SBATCH --time=00:20:00
#SBATCH --nodes=5 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/parse_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j parse_BAMs.sh

set echo on

RT="$HOME/packages/tart"

DATECODE="2023-09-29_00-47-50"
DSET="s_spcc6803"

BAM_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta500/alignments_$DATECODE"
SAVE_ROOT="$RT/validation/alignment/outputs/$DSET/plots/picks"
BOUNDS_PATH="$RT/validation/alignment/outputs/$DSET/rowid_to_bounds.json"

mpiexec tart-targeted parse_bam -i $BAM_DIR -o $SAVE_ROOT --bounds-file $BOUNDS_PATH --picks &&
    wait

echo "Finished BAM parsing for $DSET into $SAVE_ROOT"
