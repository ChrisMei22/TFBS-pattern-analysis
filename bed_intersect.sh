#!/bin/bash

# Check if the correct number of arguments is provided
# file A is TFBS from TOBIAS and file B is genomic region of interest
if [ "$#" -ne 2 ]; then
    echo "Usage: ./bed_intersect.sh <fileA> <fileB>"
    exit 1
fi

# Extract the file paths from the command-line arguments
fileA="$1"
fileB="$2"

# Define the output file path
outputFile="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Script_outputs/intersect_output.bed"

# Perform the intersection using bedtools
# See documentation on bedtools intersect to modify the options
bedtools intersect -wa -a "$fileA" -b "$fileB" | sort -u -k1,1 -k2,2 -k3,3 > "$outputFile"

# Count the number of lines in the output file
lineCount=$(wc -l < "$outputFile")
echo "Number of lines in the output file: $lineCount"
