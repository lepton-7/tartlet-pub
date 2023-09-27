#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --nodes=5 --cpus-per-task=6

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/download_SRAs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j download_SRAs.sh

# set -x
set echo on

LISTNAME="a_fischeri_ES114_griend2023_acc"
TAX="a_fischeri_ES114"

DUMP_DIR="/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/$TAX"
TEMP_DIR="/fs/scratch/PDS0325/${TAX}_temp"
ACC_LIST="/users/PDS0325/sachitk26/packages/tart/validation/rna_seq/acc_lists/$LISTNAME.txt"

module reset

module load miniconda3
module load sratoolkit/3.0.2

source activate local

echo "Downloading sequences from $(basename $ACC_LIST) into $(basename $DUMP_DIR)/"
mpiexec python -u download_SRAs.py $ACC_LIST $DUMP_DIR $TEMP_DIR &&
    wait

echo "Done"
