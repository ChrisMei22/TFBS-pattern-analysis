#!/bin/bash

# This script is to create a map of TFBS and their merged clusters with PhastCons conservation data across 15 species

# Ensure exactly two inputs are provided (TFBS and PhastCons coordinates)
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <TFBS BED file> <Conserved regions BED file>"
    exit 1
fi

# Assign command line arguments to variables
TFBS_bed="$1"
phastcons_bed="$2"

# Define paths for histone mark files
H3K4me1_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K4me1_c14a.bed"
H3K27ac_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K27ac_c14a.bed"
H3K4me3_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K4me3_c14a.bed"

# Define intermediate and output file paths
base_dir="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/cluster_stringent_histone"

# Define file paths relative to the base directory
TFBS_cluster_raw="$base_dir/TFBS_cluster_raw.bed"
TFBS_cluster="$base_dir/TFBS_cluster.bed"
merged_TFBS="$base_dir/merged_TFBS.bed"
cluster_regions="$base_dir/cluster_regions.bed"
cluster_regions_filtered="$base_dir/cluster_regions_filtered.bed"
merged_TFBS_filtered="$base_dir/merged_TFBS_filtered.bed"
merged_phastcons="$base_dir/merged_phastcons.bed"
conserved_cluster_TFBS="$base_dir/conserved_cluster_TFBS.bed"
conserved_cluster_regions="$base_dir/conserved_cluster_regions.bed"
not_conserved_cluster_TFBS="$base_dir/not_conserved_cluster_TFBS.bed"
not_conserved_cluster_regions="$base_dir/not_conserved_cluster_regions.bed"
TFBS_phastcons="$base_dir/TFBS_phastcons.bed"
cluster_phastcons="$base_dir/cluster_phastcons.bed"
complete_cluster_TFBS_conservation_map="$base_dir/complete_cluster_TFBS_conservation_map.bed"
complete_cluster_region_conservation_map="$base_dir/complete_cluster_region_conservation_map.bed"

# Step 1: Cluster TFBS according to parameters (distance = n) discovered in WBS and get only the relevant columns(18 is the one that contains clusterID):
# chrX start end clusterID
sort -k1,1 -k2,2n "$TFBS_bed" | bedtools cluster -d 12 | \
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$18}' > "$TFBS_cluster_raw"

# Step 2: get only clusters and their TFBS coordinates that have 7 or more TFBS
awk '{cluster_counts[$5]++} END {for (cluster in cluster_counts) if (cluster_counts[cluster] >=14) print cluster}' "$TFBS_cluster_raw" > "valid_clusters.txt"
while read -r cluster; do
    awk -v cluster_id="$cluster" '$5 == cluster_id' "$TFBS_cluster_raw" >> "$TFBS_cluster"
done < "valid_clusters.txt"
rm "valid_clusters.txt"

# Step 3: Merge overlapping TFBS in original clustering but keeping clusterID
sort -k1,1 -k2,2n "$TFBS_cluster" |bedtools merge -c 5 -o min > "$merged_TFBS"

# Step 4: Create the cluster regions by merging them this time using distance n+1 and keeping cluster ID info
sort -k1,1 -k2,2n "$TFBS_cluster" |bedtools merge -d 12 -c 5 -o min > "$cluster_regions"

# Step 5.1: Intersect cluster regions with enhancer histone marks (H3K4me1 AND H3K27ac)
bedtools intersect -wa -a "$cluster_regions" -b "$H3K4me1_bed" | \
bedtools intersect -wa -a - -b "$H3K27ac_bed" > "$cluster_regions_filtered"

# Step 5.2: Intersect cluster regions with promoter histone marks (H3K4me1)
#bedtools intersect -wa -a "$cluster_regions" -b "$H3K4me3_bed" >

# Step 6: Intersect merged_TFBS with filtered_cluster_regions to get non-overlapping TFBS that belong to the clusters
bedtools intersect -wa -a "$merged_TFBS" -b "$cluster_regions_filtered" > "$merged_TFBS_filtered"

# Step 7: Process PhastCons bed file - Merge overlapping conserved regions and keep the lower score of the overlap (not that it really matters)
bedtools merge -i "$phastcons_bed" -c 5 -o min > "$merged_phastcons"

bedtools intersect -a "$merged_phastcons" -b "$merged_TFBS_filtered" > "$TFBS_phastcons"
bedtools intersect -a "$merged_phastcons" -b "$cluster_regions_filtered" > "$cluster_phastcons"

# Step 8: Intersect merged conserved regions with TFBS_cluster_filtered and filtered_TFBS_regions
bedtools intersect -wa -wb -a "$merged_TFBS_filtered" -b "$TFBS_phastcons" > "$conserved_cluster_TFBS"
bedtools intersect -wa -wb -a "$cluster_regions_filtered" -b "$cluster_phastcons" > "$conserved_cluster_regions"
 
# Step 9: Get the regions and TFBS that do not overlap with PhastCons and concatenate both files to get complete map:
bedtools intersect -v -a "$merged_TFBS_filtered" -b "$TFBS_phastcons" > "$not_conserved_cluster_TFBS"
bedtools intersect -v -a "$cluster_regions_filtered" -b "$cluster_phastcons" > "$not_conserved_cluster_regions"

cat "$conserved_cluster_TFBS" "$not_conserved_cluster_TFBS" > "$complete_cluster_TFBS_conservation_map"
cat "$conserved_cluster_regions" "$not_conserved_cluster_regions" > "$complete_cluster_region_conservation_map"

echo "Pipeline completed. Outputs are located at:"
echo "Output 1: $complete_cluster_TFBS_conservation_map"
echo "Output 2: $complete_cluster_region_conservation_map"

