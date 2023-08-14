#!/bin/bash

#SBATCH --time=01:30:00
#SBATCH --nodes=1 --cpus-per-task=44

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# There are 92 pairs of readsets
#SBATCH --array=0-92

# RUN THIS TO LAUNCH JOB
## DATECODE="$(date +"%Y-%m-%d_%H-%M-%S")"; sbatch --output=jobs/align_reads.out.$DATECODE.%j align_reads.sh

# set -x
set echo on

module load hisat2

DATECODE="2023-08-14_11-02-43"

IDX=$((SLURM_ARRAY_TASK_ID))

RT="$HOME/packages/tart"

REF_DIR="$RT/validation/alignment/switch_seqs_delta500"
SEQ_DIR="$RT/validation/rna_seq/e_coli"

mkdir -p $REF_DIR

ALIGN_DIR="$REF_DIR/alignments_$DATECODE"
mkdir -p $ALIGN_DIR

# Look for all the paired sequence data
cd $SEQ_DIR
FQ_FILES=(*_1.fastq.gz)

# Pick one pair to align with per script
SEQ_FILE="$SEQ_DIR/${FQ_FILES[$IDX]:0:-11}"

cd $REF_DIR

REFS=(*.fna)

# Iterate over each riboswitch class reference file
for Q in ${REFS[@]}; do

    RSWTCH_CLASS=${Q:0:-4}

    IDX_DIR="${RSWTCH_CLASS}_index"

    # Directory for SAM output
    SAMOUT_DIR=$ALIGN_DIR/$RSWTCH_CLASS
    mkdir -p $SAMOUT_DIR

    hisat2 -p 44 -x $IDX_DIR/$IDX_DIR \
        -1 "$SEQ_FILE"_1.fastq.gz -2 "$SEQ_FILE"_2.fastq.gz \
        -S $SAMOUT_DIR/$RSWTCH_CLASS.$(basename "$SEQ_FILE").sam -t \
        --no-unal --score-min L,0,-0.6

    echo
    echo "---------------------------------------------------------------------"
    echo
done
