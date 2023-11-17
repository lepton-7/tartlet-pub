#!/bin/bash

#SBATCH --time=03:00:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# --------------------------------
# There are 24 pairs of readsets
#SBATCH --array=0-24
# --------------------------------

# RUN THIS TO LAUNCH JOB
## DATECODE="$(date +"%Y-%m-%d_%H-%M-%S")"; sbatch --output=jobs/4.m_tuber.align_reads.out.$DATECODE.%j 4.align_reads.sh

# set -x
set echo on

IDX=$((SLURM_ARRAY_TASK_ID))

DATECODE="20231113"
DSET="m_tuber"

RT="$HOME/packages/tart"

REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000"
SEQ_DIR="$RT/validation/rna_seq/$DSET"

mkdir -p $REF_DIR

ALIGN_DIR="$REF_DIR/alignment_$DATECODE"
mkdir -p $ALIGN_DIR

# Look for all the paired sequence data
cd $SEQ_DIR
FQ_FILES=(*_1.fastq.gz)

# Pick one pair to align with per script
SEQ_FILE="$SEQ_DIR/${FQ_FILES[$IDX]:0:-11}"

cd $REF_DIR

echo "Instance $IDX starting alignment for $DSET from $SEQ_DIR reads"
echo
echo

srun tart-targeted align -i $REF_DIR \
    -o $ALIGN_DIR \
    -1 "$SEQ_FILE"_1.fastq.gz \
    -2 "$SEQ_FILE"_2.fastq.gz \
    --readpair-name $(basename "$SEQ_FILE") \
    -t \
    --no-unal \
    --score-min L,0,-0.6

echo "Finished $(basename "$SEQ_FILE") read alignment for $DSET into $ALIGN_DIR"
