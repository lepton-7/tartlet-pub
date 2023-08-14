#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/make_index.out.$(date +"%Y-%m-%d_%H-%M-%S").%j make_index.sh

# set -x
set echo on

module load hisat2

RT="$HOME/packages/tart"
REF_DIR="$RT/validation/alignment/switch_seqs_delta500"

cd $REF_DIR

arr=(*.fna)

for i in ${arr[@]}; do

    iPATH="${i:0:-4}_index"

    mkdir -p $iPATH

    hisat2-build -p 44 $i $iPATH/${i:0:-4}_index
done
