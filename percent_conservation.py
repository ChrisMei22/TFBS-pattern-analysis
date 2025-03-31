# This script will create the scatterplot to map percentage conservation for my candidate clusters
import matplotlib.pyplot as plt


# Files in this script
tfbs_cluster_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_TFBS_conservation_map.bed"
region_cluster_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_region_conservation_map.bed"

# Process files into a dictionary in the following format -> clusterID: [total bp, conserved bp]
def extract_TFBS(file_path: str) -> dict:
    cluster_dict = {}  # Initialize dictionary to hold clusterID: [total bp, conserved bp]

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()

            # Extract necessary columns, default to 0 for conserved bp if columns 5 and 6 are missing
            total_bp = int(split_line[2]) - int(split_line[1])
            conserved_bp = 0  # Default value

            # Check if conserved region columns exist
            if len(split_line) >= 8:
                conserved_bp = int(split_line[6]) - int(split_line[5])

            cluster_id = split_line[3]

            # Update the dictionary
            if cluster_id in cluster_dict:
                cluster_dict[cluster_id][0] += total_bp
                cluster_dict[cluster_id][1] += conserved_bp
            else:
                cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict

def merge_intervals(intervals):
    """Merge overlapping intervals."""
    if not intervals:
        return []

    # Sort intervals based on start position
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]

    for current_start, current_end in intervals[1:]:
        last_end = merged[-1][1]

        # If the current interval overlaps with the last merged interval, use the max end value
        if current_start <= last_end:
            merged[-1] = (merged[-1][0], max(last_end, current_end))
        else:
            merged.append((current_start, current_end))
    return merged

def calculate_total_bp(intervals):
    """Calculate total base pairs from merged intervals."""
    return sum(end - start for start, end in intervals)

def extract_region(file_path: str) -> dict:
    cluster_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            cluster_id = split_line[3]
            total_interval = (int(split_line[1]), int(split_line[2]))
            conserved_interval = (int(split_line[5]), int(split_line[6])) if len(split_line) >= 7 else None

            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = {'total': [total_interval], 'conserved': [] if conserved_interval is None else [conserved_interval]}
            else:
                cluster_dict[cluster_id]['total'].append(total_interval)
                if conserved_interval:
                    cluster_dict[cluster_id]['conserved'].append(conserved_interval)

    # Merge intervals and calculate bp
    for cluster_id, intervals in cluster_dict.items():
        merged_total = merge_intervals(intervals['total'])
        merged_conserved = merge_intervals(intervals['conserved'])

        total_bp = calculate_total_bp(merged_total)
        conserved_bp = calculate_total_bp(merged_conserved)

        cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict


tfbs_dict = extract_TFBS(tfbs_cluster_file)
#print(tfbs_dict)
region_dict = extract_region(region_cluster_file)
#print(region_dict)
tfbs_percentage_array = []
region_percentage_array = []
key_array = []

conservation_map = {}
tfbs_percent_conserved = 0
region_percent_conserved = 0
for key, values in tfbs_dict.items():
    tfbs_percent_conserved = values[1]/values[0] * 100
    tfbs_percentage_array.append(tfbs_percent_conserved)
    key_array.append(key)
    if key in region_dict.keys():
        region_percent_conserved = region_dict[key][1]/region_dict[key][0] * 100
        region_percentage_array.append(region_percent_conserved)
    conservation_map[key] = (tfbs_percent_conserved, region_percent_conserved)

print("total number of clusters: ", len(key_array))
print (conservation_map)

# THis is to find any mismatches
tfbs_keys = set(tfbs_dict.keys())
region_keys = set(region_dict.keys())

# Find keys that are in tfbs_dict but not in region_dict
unmatched_tfbs_keys = tfbs_keys - region_keys

# Find keys that are in region_dict but not in tfbs_dict
unmatched_region_keys = region_keys - tfbs_keys

print("Unmatched TFBS keys:", unmatched_tfbs_keys)
print("Unmatched Region keys:", unmatched_region_keys)

plt.scatter(region_percentage_array,tfbs_percentage_array, c="#e69138")
plt.xlabel("% cluster region bp conserved")
plt.ylabel("% TFBS bp conserved")
plt.title("cCRE conservation")
plt.show()

# Histogram with total region bp conservation since relationship is linear
plt.hist(region_percentage_array, color="#e69138", edgecolor = "white")
plt.xlabel("% conservation")
plt.ylabel("Frequency")
plt.title("cCRE conservation")
plt.axvline(70, color='k', linestyle='dashed', linewidth=1)
plt.show()

high_conservation_cluster: list = []
for key,values in conservation_map.items():
    if values[1] >= 70:
        high_conservation_cluster.append(key)
print(high_conservation_cluster)
