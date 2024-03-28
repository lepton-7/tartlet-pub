DSETS=(
    "a_fischeri_ES114"
    "a_kunk"
    # "b_anth"
    # "b_frag"
    # "b_pseudo"
    # "b_sub_168"
    # "b_theta"
    # "b_xyla"
    # "c_diff"
    # "c_vibrioides"
    # "d_vulg"
    # "e_coli"
    # "e_fae"
    # "e_limo"
    # "k_pneum"
    # "m_smeg"
    # "m_tuber"
    # "n_gonorr"
    # "p_cholor_aureo3084"
    # "p_fluor"
    # "p_salmo"
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
# These were removed from ^
#   "c_basil"
#   "p_aeru"

for DSET in ${DSETS[@]}; do
    DATECODE="$(date +"%Y-%m-%d_%H-%M-%S")"
    SEQ_DIR="/fs/ess/PDS0325/Riboswitches/data/rna_seq/$DSET"
    
    PAIRS=$(ls $SEQ_DIR/*_1.fastq.gz | wc -l)

    sbatch --export=ALL,DSET="$DSET" --array=0-$PAIRS --output=./4.$DSET.align_reads.out.$DATECODE.%j test.sh
done
