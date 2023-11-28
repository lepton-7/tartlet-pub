#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=2 --cpus-per-task=2

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# Change the output file to reflect input directory used
## sbatch --output=jobs/1.call_orfs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 1.call_orfs.sh

set echo on

RT="$HOME/packages/tart"
GENOMES="$RT/validation/genomes"
OUT="$RT/validation/orfs/calls"

mpiexec tart-utils find-orfs -o $OUT $GENOMES/*.fna &&
    wait

echo "Done"
