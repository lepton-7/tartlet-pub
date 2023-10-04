#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/find_switch_seqs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j find_switch_seqs.sh

#
#
#
# Target script reworked to targeted/make_reference_sequences.py; REFACTOR SCRIPT
#
#
#

set -x
set echo on

module load miniconda3
source activate local

RT="$HOME/packages/tart"
TABLE="$RT/validation/tables/inf_results.csv"
GENOMES="$RT/validation/genomes"

DELTA=500
DSET="s_spcc6803"

OUT_DIR="$RT/validation/alignment/outputs/$DSET"

mpiexec python -u find_switch_seqs.py $TABLE $GENOMES $DELTA $DSET $OUT_DIR &&
    wait
