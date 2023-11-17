#!/bin/bash

#SBATCH --time=00:50:00
#SBATCH --nodes=2 --cpus-per-task=1

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# RUN THIS TO LAUNCH JOB
## sbatch --output=jobs/6.all.parse_BAMs.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 6.parse_BAMs.sh

set echo on

RT="$HOME/packages/tart"

DATECODE="20231113"

DSETS=(
    "a_baum"
    "a_rabiei"
    "c_vibrioides"
    "m_tuber"
    "s_elon"
    "s_spcc6803"
    "a_fischeri_ES114"
    "b_sub_168"
    "e_coli"
    "p_fluor"
    "s_enter_typh"
    "x_albi"
    "a_kunk"
    "b_theta"
    "k_pneum"
    "s_coelicolor"
    "s_sanguinis"
)

for DSET in ${DSETS[@]}; do

    BAM_DIR="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000/alignment_$DATECODE"
    SAVE_ROOT="$RT/validation/alignment/outputs/$DSET/plots/picks"
    BOUNDS_PATH="$RT/validation/alignment/outputs/$DSET/switch_seqs_delta1000-1000/rowid_to_bounds.json"

    echo "Starting BAM parsing for $DSET"

    mpiexec tart-targeted parse-bam -i $BAM_DIR -o $SAVE_ROOT --bounds-file $BOUNDS_PATH --picks --allow-soft-clips &&
        wait

    echo "Finished BAM parsing for $DSET into $SAVE_ROOT"
    echo "----------------------------------------------"
    echo
done
