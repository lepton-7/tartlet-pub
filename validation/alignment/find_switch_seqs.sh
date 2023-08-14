#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/find_switch_seqs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j find_switch_seqs.sh

set -x
set echo on

module load miniconda3
source activate local

mpiexec python -u find_switch_seqs.py ../complete_tax_downstream.csv $HOME/Emerge_MAGs_v1/Derep95 &&
    wait
