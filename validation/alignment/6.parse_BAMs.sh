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
    # "a_baum"
    # "a_fischeri_ES114"
    # "a_kunk"
    # "a_rabiei"
    # "b_anth"
    # "b_frag"
    # "b_pseudo"
    "b_sub_168"
    # "b_theta"
    # "b_xyla"
    # "c_basil"
    # "c_diff"
    # "c_vibrioides"
    # "d_vulg"
    "e_coli"
    # "e_fae"
    # "e_limo"
    # "k_pneum"
    # "m_smeg"
    # "m_tuber"
    # "n_gonorr"
    # "p_fluor"
    # "s_coelicolor"
    # "s_elon"
    # "s_enter_typh"
    # "s_epi"
    # "s_meli"
    # "s_sanguinis"
    # "s_spcc6803"
    # "x_albi"
    # "x_ory"
)

echo "Parsing BAMs for:"
echo "${DSETS[*]}"
echo
echo

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

echo "Done"
