#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --nodes=1 --cpus-per-task=40

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## bash ./4_1.launch_alignarray.sh

# DEPRECATED
## DATECODE="$(date +"%Y-%m-%d_%H-%M-%S")"; sbatch --output=jobs/4.b_sub_168.align_reads.out.$DATECODE.%j 4.align_reads.sh

# set -x
set echo on

IDX=$((SLURM_ARRAY_TASK_ID))

DATECODE="20240328"
# DSET="p_salmo"

RT="$HOME/packages/tart"

REF_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000"
SEQ_DIR="/fs/ess/PDS0325/Riboswitches/data/rna_seq/$DSET"

mkdir -p $REF_DIR

ALIGN_DIR="$REF_DIR/alignment_$DATECODE"
mkdir -p $ALIGN_DIR

# Look for all the paired sequence data
cd $SEQ_DIR
FQ_FILES=(*_1.fastq.gz)

# Pick one pair to align with per script; has extension .fastq.gz
# But want to remove _n.fastq.gz
SEQ_FILE="$SEQ_DIR/${FQ_FILES[$IDX]:0:-11}"

cd $REF_DIR

echo "Instance $IDX starting alignment for $DSET from $SEQ_DIR reads"
echo
echo

# Only want to run on the unified set
FSTA="unified.fna"
STEM=${FSTA:0:-4}

mkdir -p "$ALIGN_DIR/$STEM"

srun hisat2 \
    -x "$REF_DIR/${STEM}_index/${STEM}_index" \
    -1 "$SEQ_FILE"_1.fastq.gz \
    -2 "$SEQ_FILE"_2.fastq.gz \
    -S "$ALIGN_DIR/$STEM/$(basename "$SEQ_FILE").sam" \
    -t \
    --no-unal \
    --score-min L,0,-0.4 \
    -p 40

# srun tart-targeted align -i $REF_DIR \
#     -o $ALIGN_DIR \
#     -1 "$SEQ_FILE"_1.fastq.gz \
#     -2 "$SEQ_FILE"_2.fastq.gz \
#     --readpair-name $(basename "$SEQ_FILE") \
#     -t \
#     --no-unal \
#     --score-min L,0,-0.4 \
#     -p 40

echo "Finished $(basename "$SEQ_FILE") read alignment for $DSET into $ALIGN_DIR"
