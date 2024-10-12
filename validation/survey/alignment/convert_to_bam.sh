DIR="/fs/scratch/PDS0325/e_coli_alignment_20241011"

cd $DIR

arr=($DIR/*.sam)

for SAM in ${arr[@]}; do
    stem=$(basename $SAM .sam)

    echo "Processing $stem"

    samtools view -@ 46 -Sb -o "$stem.bam" "$stem.sam"
    samtools sort -@ 46 -O bam -o "$stem.sorted.bam" "$stem.bam"
    samtools index "$stem.sorted.bam"

    echo "Finished processing $stem"
    echo
done
