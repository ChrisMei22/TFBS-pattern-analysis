#!/bin/bash
#mounts USB disk 
diskutil mount /Volumes/MEI_DROS

cd /Volumes/MEI_DROS/Drosophila_Materials_Drive/Drosophila_Dataset

# template: for dir in ~/projects/git/*; do (cd "$dir" && git pull); done
# The code inside the bigger loop seems to work, the problem is that the script doesn't go into the folders

for dir in */ 
do
    (echo $dir
    cd $dir
    for file_gff in *.gff3
    do
        gff2bed < $file_gff > "$(basename $file_gff .gff3)".bed
        sed -e 's/^Y/chrY/' -e 's/^X/chrX/' -e 's/^2L/chr2L/' -e 's/^2R/chr2R/'  -e 's/^3L/chr3L/' -e 's/^3R/chr3R/' -e 's/^4/chr4/' \
        -e 's/^M/chrM/' "$(basename $file_gff .gff3)".bed > "$(basename $file_gff .gff3)".chr.bed
    done
    # fix bed and bam files to chr versions + delete the bam.bai file
    for file_bam in *.bam; do
        samtools view -H $file_bam | \
                sed -e 's/SN:2L/SN:chr2L/' | sed -e 's/SN:2R/SN:chr2R/' | \
                sed -e 's/SN:3L/SN:chr3L/' | sed -e 's/SN:3R/SN:chr3R/' | \
                sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
                sed -e 's/SN:M/SN:chrM/' | samtools reheader - $file_bam > ./"$(basename $file_bam)".chr.bam
    #if loop inside for loop to remove bam.bai file
        rm *.bai
    done)
done
