#!/bin/bash

for dir in */; do
    cd "$dir" &&
    find . -type f -name "*.bai" -delete
    TOBIAS ATACorrect --bam "$(find . -name "*_chr.bam*")" --genome /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep1/final_dm3.fa --peaks "$(find . -name "*chr.bed")"  --blacklist /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep2/dm3-blacklist.v2.bed --outdir ./ATACorrect_test --cores 8
    
    TOBIAS FootprintScores --signal "$(find . -name "*_corrected.bw*")" --regions "$(find . -name "*chr.bed")" --output ./footprints.bw --cores 8
    
    TOBIAS BINDetect --motifs /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep2/motif2.meme --motif_pvalue 0.001 --signals footprints.bw --genome /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep1/final_dm3.fa --peaks "$(find . -name "*chr.bed")"  --outdir ./BINDetect_output --cond_names "$(basename "$dir")" --cores 8
    cd -
done



