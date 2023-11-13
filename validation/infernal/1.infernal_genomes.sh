#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=2 --cpus-per-task=2

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# Change the output file to reflect input directory used
## sbatch --output=jobs/1.infernal_genomes.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 1.infernal_genomes.sh

set echo on

RT="$HOME/packages/tart"
GENOMES="$RT/validation/genomes"
OUT="$RT/validation/infernal/results"

mpiexec tart-utils find-riboswitches -o $OUT --no-stats $GENOMES/*.fna &&
    wait

echo "Done"
