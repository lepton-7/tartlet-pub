#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/parse_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 6.parse_BAMs.sh

set echo on

RT="$HOME/packages/tart"

DATECODE="2023-08-14_11-02-43"
DSET="e_coli"

BAM_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta500/alignments_$DATECODE"
SAVE_ROOT="$RT/validation/alignment/outputs/$DSET/plots/picks"
BOUNDS_PATH="$RT/validation/alignment/outputs/$DSET/rowid_to_bounds.json"

mpiexec tart-targeted parse-bam -i $BAM_DIR -o $SAVE_ROOT --bounds-file $BOUNDS_PATH --picks --allow-soft-clips &&
    wait

echo "Finished BAM parsing for $DSET into $SAVE_ROOT"
