#!/bin/bash

# Define the directory to search and the output directory
input_dir="/users/omokhtar/PDBbind" # folder of subfolders
output_dir="/users/omokhtar/representatives"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through the directories in the input directory
for folder in "$input_dir"/*; do
  if [ -d "$folder" ]; then
    # Extract the folder name
    folder_name=$(basename "$folder")
    
    # Define the source and destination paths
    src_file="$folder/analysis/clusters.pdb"
    dest_file="$output_dir/${folder_name}.pdb"
    
    # Check if the source file exists, then copy it
    if [ -f "$src_file" ]; then
      cp "$src_file" "$dest_file"
    fi
  fi
done

