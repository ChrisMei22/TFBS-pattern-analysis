# This script is to measure the % conservation for the region as a whole for my 5 candidate enhancers
# File for this study is a map based on bedtools intersects between cCRM regions and overlapping Phastcons chunks

ccrm_conservation_map_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/candidate_enhancer_conservation_map.txt"
ccrm_conservation: dict = {}
with open(ccrm_conservation_map_file, 'r') as file:
    for line in file:
        split_line = line.strip().split()
        conserved_bp = 0

        if len(split_line) >= 9:
            conserved_bp = int(split_line[6]) - int(split_line[5])

        cCRM_ID = split_line[3]
        cCRM_size : int = int(split_line[2]) - int(split_line[1])
        cCRM_key = (cCRM_ID, cCRM_size)
        if cCRM_key in ccrm_conservation:
            ccrm_conservation[cCRM_key] += conserved_bp
        else:
            ccrm_conservation[cCRM_key] = conserved_bp
print(ccrm_conservation)

for key, value in ccrm_conservation.items():
    percentage = (value / key[1])*100
    print("Cluster:", key[0], "has a conservation of", percentage, "%")




