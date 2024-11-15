#!/bin/bash

#SBATCH --time=00:20:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/fragment_lengths.out.$(date +"%Y-%m-%d_%H-%M-%S").%j fragment_lengths.sh

set echo on

SCR="/fs/scratch/PDS0325"
RT="$HOME/tartlet-pub"

INPUT_DIR="/fs/scratch/PDS0325/e_coli_alignment_20241011"
OUTPUT_PATH="e_coli_fragment_sizes_2.csv"

echo "Initialising."

mpiexec python fragment_lengths.py -i $INPUT_DIR -o $OUTPUT_PATH &&
    wait

echo
echo
echo "Done"
