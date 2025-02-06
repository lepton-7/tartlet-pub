#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=6 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/cophen_d_variations.out.$(date +"%Y-%m-%d_%H-%M-%S").%j cophen_d_variations.sh

set echo on

mpiexec python -u cophen_d_variations.py

echo "Done"
