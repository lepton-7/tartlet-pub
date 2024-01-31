#!/bin/bash

#SBATCH --time=00:20:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/fragment_lengths.out.$(date +"%Y-%m-%d_%H-%M-%S").%j fragment_lengths.sh

set echo on

SCR="/fs/scratch/PDS0325"
RT="$HOME/packages/tart"

INPUT_DIR="$RT/validation/survey/alignment/b_sub_168_alignment_20231113"
OUTPUT_PATH="b_sub_168_fragment_sizes.csv"

echo "Initialising."

mpiexec python fragment_lengths.py -i $INPUT_DIR -o $OUTPUT_PATH &&
    wait

echo
echo
echo "Done"
