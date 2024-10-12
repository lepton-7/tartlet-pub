#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# --------------------------------
# There are 92 pairs of readsets
#SBATCH --array=0-92
# --------------------------------

# RUN THIS TO LAUNCH JOB
## DATECODE="$(date +"%Y-%m-%d_%H-%M-%S")"; sbatch --output=../jobs/e_coli.align.out.$DATECODE.%j align.sh

set echo on

IDX=$((SLURM_ARRAY_TASK_ID))

DATECODE="20241011"
DSET="e_coli"

RT="$HOME/tartlet-pub"

RUN_DIR="$RT/validation/survey/alignment"
IDX_DIR="$RT/validation/survey/${DSET}_index"
SEQ_DIR="/fs/ess/PDS0325/Riboswitches/data/rna_seq/$DSET"

# ALIGN_DIR="$RUN_DIR/${DSET}_alignment_$DATECODE"
ALIGN_DIR="/fs/scratch/PDS0325/${DSET}_alignment_$DATECODE"

mkdir -p $ALIGN_DIR

# Look for all the paired sequence data
cd $SEQ_DIR
FQ_FILES=(*_1.fastq.gz)

# Pick one pair to align with per script
SEQ_FILE="$SEQ_DIR/${FQ_FILES[$IDX]:0:-11}"

cd $RUN_DIR

echo "Instance $IDX starting alignment for $DSET from $SEQ_DIR reads"
echo
echo

srun hisat2 \
    -x $IDX_DIR/${DSET}_index \
    -p 44 \
    -S $ALIGN_DIR/$(basename "$SEQ_FILE").sam \
    -1 "$SEQ_FILE"_1.fastq.gz \
    -2 "$SEQ_FILE"_2.fastq.gz \
    -t \
    --no-unal \
    --score-min L,0,-0.6

echo "Finished $(basename "$SEQ_FILE") read alignment for $DSET into $ALIGN_DIR"
