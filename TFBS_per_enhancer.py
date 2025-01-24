# This script is to count the number of TFBS per known enhancers
import matplotlib.pyplot as plt
import numpy as np

def tfbs_per_CRM (file_path: str) -> dict:
    tfbs_crm_map = {}
    no_dups_list = []
    with open(file_path, 'r') as file:
        raw_list = []
        for line in file:
            split_line = line.strip().split()
            processed_line = []
            for index, item in enumerate(split_line):
                if index in [1, 2, 5, 6]:
                    processed_line.append(int(item))  # Convert to integer
                else:
                    processed_line.append(item)  # Keep as string
            raw_list.append(processed_line)
    # Removing duplicates
    for item in raw_list:
        if item not in no_dups_list:
            no_dups_list.append(item)
    print("Total unique combinations: ", len(no_dups_list))

    for item in no_dups_list:
        crm_id: tuple = (item[1], item[2], item[3])
        if crm_id not in tfbs_crm_map:
            tfbs_crm_map[crm_id] = 1
        else:
            tfbs_crm_map[crm_id] += 1

    return tfbs_crm_map

def tfbs_count_dist (tfbs_CRM_map: dict) -> None:
    tfbs_count = []
    for tf_count in tfbs_CRM_map.values():
        tfbs_count.append(tf_count)
    print ("Mean: ", np.mean(tfbs_count))
    print ("Median: ", np.median(tfbs_count))

    plt.hist(tfbs_count,bins = 20,color = '#db5856', edgecolor = 'white')
    plt.xlabel("TFBS per known CRM")
    plt.ylabel("Frequency")
    plt.title("TFBS distribution per known CRM (stage 4-6)")
    plt.axvline(np.median(tfbs_count), color='k', linestyle='dashed', linewidth=1)
    plt.show()
    return None

tfbs_crm_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/CRM_TFBS_map.txt"

tfbs_count_dist(tfbs_per_CRM(tfbs_crm_file))
