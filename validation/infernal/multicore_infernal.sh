#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --nodes=1 --cpus-per-task=2

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# Change the output file to reflect input directory used
## sbatch --output=jobs/multicore_infernal.out.$(date +"%Y-%m-%d_%H-%M-%S").%j multicore_infernal.sh

set -x
set echo on

module use $HOME/local/share/lmodfiles
module load infernal
module load miniconda3

source activate local

RT="$HOME/packages/tart"
GENOMES="$RT/validation/genomes"
OUT="$RT/validation/infernal/results"

rswitches="$HOME/cm/riboswitches/riboswitches.cm"
rfam="$HOME/cm/Rfam/Rfam.cm"

mpiexec python multicore_infernal.py $GENOMES $OUT $rswitches &&
    wait
