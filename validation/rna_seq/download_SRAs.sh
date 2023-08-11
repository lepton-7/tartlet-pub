#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/download_SRAs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j download_SRAs.sh

set -x
set echo on

DUMP_DIR="/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/e_coli"
ACC_LIST="/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/e_coli_sastry2019_acc.txt"

module load miniconda3
source activate local

module load sratoolkit/3.0.2

mpiexec python -u download_SRAs.py $ACC_LIST $DUMP_DIR


