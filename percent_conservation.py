"""
Script Name: percent_conservation.py
Description: Compares and visualizes the evolutionary conservation profiles of individual 
             TFBS tracks versus their encompassing macro-cluster regions. Generates 
             scatter plots, distribution histograms, and filters
             high-conservation candidate clusters that have passed the clustering and 
             histone mark criteria based on a 70% overlap threshold.
Author: Christian Mei
"""

import matplotlib.pyplot as plt

# --- GLOBAL FILE PATH DEFINITIONS ---
tfbs_cluster_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_TFBS_conservation_map.bed"
region_cluster_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_region_conservation_map.bed"


# --- DATA PARSING & PROCESSING FUNCTIONS ---

def extract_TFBS(file_path: str) -> dict:
    """
    Parses individual TFBS conservation files to map absolute footprint lengths 
    and matching PhastCons base pair counts to their respective cluster IDs.

    Parameters:
        file_path (str): Path to the TFBS conservation BED map.

    Returns:
        dict: Format -> clusterID: [total_bp, conserved_bp]
    """
    cluster_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()

            # Calculate the absolute window length of the individual TFBS feature
            total_bp = int(split_line[2]) - int(split_line[1])
            conserved_bp = 0  

            # Accumulate conservation width if a PhastCons intersection track exists
            if len(split_line) >= 8:
                conserved_bp = int(split_line[6]) - int(split_line[5])

            cluster_id = split_line[3]

            # Aggregate overlapping footprint entries under a unified cluster identifier
            if cluster_id in cluster_dict:
                cluster_dict[cluster_id][0] += total_bp
                cluster_dict[cluster_id][1] += conserved_bp
            else:
                cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict


def merge_intervals(intervals: list) -> list:
    """
    Collapses overlapping or redundant coordinate intervals to ensure precise 
    base pair quantification without multi-counting intersection boundaries.

    Parameters:
        intervals (list): Unsorted list of (start, end) coordinate tuples.

    Returns:
        list: Consolidated list of sorted, non-overlapping coordinate tuples.
    """
    if not intervals:
        return []

    # Sort intervals hierarchically based on structural start positions
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]

    for current_start, current_end in intervals[1:]:
        last_end = merged[-1][1]

        # Resolve overlaps by extending the terminal boundary of the previous block
        if current_start <= last_end:
            merged[-1] = (merged[-1][0], max(last_end, current_end))
        else:
            merged.append((current_start, current_end))
    return merged


def calculate_total_bp(intervals: list) -> int:
    """Calculates cumulative nucleotide length across processed coordinate intervals."""
    return sum(end - start for start, end in intervals)


def extract_region(file_path: str) -> dict:
    """
    Parses macro-cluster tracking files, applying interval deduplication 
    to calculate precise global conservation ratios across entire loci.

    Parameters:
        file_path (str): Path to the macro-cluster region conservation BED map.

    Returns:
        dict: Format -> clusterID: [total_region_bp, conserved_region_bp]
    """
    cluster_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            cluster_id = split_line[3]
            total_interval = (int(split_line[1]), int(split_line[2]))
            conserved_interval = (int(split_line[5]), int(split_line[6])) if len(split_line) >= 7 else None

            # Track raw spatial intervals separately prior to deduplication
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = {
                    'total': [total_interval], 
                    'conserved': [] if conserved_interval is None else [conserved_interval]
                }
            else:
                cluster_dict[cluster_id]['total'].append(total_interval)
                if conserved_interval:
                    cluster_dict[cluster_id]['conserved'].append(conserved_interval)

    # Resolve overlapping regions and evaluate absolute nucleotide widths
    for cluster_id, intervals in cluster_dict.items():
        merged_total = merge_intervals(intervals['total'])
        merged_conserved = merge_intervals(intervals['conserved'])

        total_bp = calculate_total_bp(merged_total)
        conserved_bp = calculate_total_bp(merged_conserved)

        cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict


# --- DATA PIPELINE & CALCULATIONS ---

tfbs_dict = extract_TFBS(tfbs_cluster_file)
region_dict = extract_region(region_cluster_file)

tfbs_percentage_array = []
region_percentage_array = []
key_array = []
conservation_map = {}

# Compute percentage ratios between regional features and footprint components
for key, values in tfbs_dict.items():
    tfbs_percent_conserved = (values[1] / values[0]) * 100
    tfbs_percentage_array.append(tfbs_percent_conserved)
    key_array.append(key)
    
    if key in region_dict:
        region_percent_conserved = (region_dict[key][1] / region_dict[key][0]) * 100
        region_percentage_array.append(region_percent_conserved)
    
    conservation_map[key] = (tfbs_percent_conserved, region_percent_conserved)

print("Total number of successfully mapped clusters:", len(key_array))
print("Conservation mapping directory entries:", conservation_map)


# --- QUALITY CONTROL MAPPING INTEGRITY CHECK ---

tfbs_keys = set(tfbs_dict.keys())
region_keys = set(region_dict.keys())

unmatched_tfbs_keys = tfbs_keys - region_keys
unmatched_region_keys = region_keys - tfbs_keys

print("Unmatched TFBS keys:", unmatched_tfbs_keys)
print("Unmatched Region keys:", unmatched_region_keys)


# --- DATA VISUALIZATION CONSTRUCTION ---

# Figure 1: Covariation Scatter Plot (TFBS vs Macro-Cluster Region Conservation)
plt.figure(figsize=(7, 6))
plt.scatter(region_percentage_array, tfbs_percentage_array, c="#e69138", alpha=0.7, edgecolors='none', s=40)
plt.xlabel("% Cluster Region bp Conserved", fontsize=11, labelpad=8)
plt.ylabel("% TFBS bp Conserved",
