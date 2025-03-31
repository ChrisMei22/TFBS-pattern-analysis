#!/bin/bash

# Find folders with 'negative_control_info' in the current directory
for folder in *_Rep1; do
    if [[ -d "$folder" ]]; then
        # Check if BINDetect_output folder exists
        bindetect_folder="$folder/BINDetect_output"
        if [[ -d "$bindetect_folder" ]]; then
            # Create an array to store the found bound.bed files
            bound_files=()
            
            # Find bound.bed files (excluding unbound.bed) in BINDetect_output folder
            while IFS= read -r -d '' file; do
                if [[ "$file" != *"unbound.bed" ]]; then
                    bound_files+=("$file")
                fi
            done < <(find "$bindetect_folder" -name "*bound.bed" -type f -print0)
            
            # Compile the bound.bed files into a single file
            output_file="${folder}_bound.bed"
            cat "${bound_files[@]}" > "$output_file"
            
            echo "Compiled file created: $output_file"
        fi
    fi
done
