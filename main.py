# This code helps turn RAMPAGE data into promoters that include 300 bp and 60bp upstream and downstream from the 5' end
rampage_file: str = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_TSC.bed"
neg_tss: list = []
pos_tss: list = []

r5_57_gff_tss_file: str = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/dmel_r5.57_gff_based_genome_features/polished_genomic_features_no_chrU/dmel-all-r5.57-TSS_chr.bed"

def read_file (input_file: str) -> list:
    txt_list = []
    with open(input_file, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            # Convert split_line[1] and split_line[2] to integers
            split_line[1] = int(split_line[1])
            split_line[2] = int(split_line[2])
            txt_list.append(split_line)
    return txt_list

# print(read_utr_file(utr_file))
rampage_file_list = read_file(rampage_file)
r5_57_gff_tss_list = read_file(r5_57_gff_tss_file)

def tss_to_promoter_300_60 (tss_list: list) -> list:
    promoter_list: list = []

    for item in tss_list:
        # Since there is a range involved, take the center of the RAMPAGE/TSS peak
        promoter_center = int((item[1] + item[2])/2)
        # Add bp according to strand direction
        if item[5] == "+":
            promoter_start = promoter_center - 300
            promoter_end = promoter_center + 60
            promoter_list.append([item[0],promoter_start, promoter_end, item[3], item[5]])
        elif item[5] == "-":
            promoter_start = promoter_center - 60
            promoter_end = promoter_center + 300
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
    return promoter_list

# Different function for promoter range defined in
def tss_to_promoter_250_50 (tss_list: list) -> list:
    promoter_list: list = []

    for item in tss_list:
        promoter_center = int((item[1] + item[2])/2)
        if item[5] == "+":
            promoter_start = promoter_center - 250
            promoter_end = promoter_center + 50
            promoter_list.append([item[0],promoter_start, promoter_end, item[3], item[5]])
        elif item[5] == "-":
            promoter_start = promoter_center - 50
            promoter_end = promoter_center + 250
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
    return promoter_list


def export_to_bed(data_list, output_file):
    with open(output_file, 'w') as file:
        for entry in data_list:
            # Join each element of the entry with a tab and write to file
            line = "\t".join(map(str, entry))
            file.write(line + "\n")
    return None

# Testing RAMPAGE data to -300 and +60 definition of promoter
#export_to_bed(tss_to_promoter_300_60(rampage_file_list), "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_promoter_300_60.bed")

# Testing RAMPAGE data to -250 and +50 definition of promoter
export_to_bed(tss_to_promoter_250_50(rampage_file_list), "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_250_50_promoter.bed")
# Testing GFF derived TSS data to -250 and +50 definition of promoter
export_to_bed(tss_to_promoter_250_50(r5_57_gff_tss_list), "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/dmel_r5.57_gff_based_genome_features/polished_genomic_features_no_chrU/dmel-all-r5.57-TSS_promoter.bed")