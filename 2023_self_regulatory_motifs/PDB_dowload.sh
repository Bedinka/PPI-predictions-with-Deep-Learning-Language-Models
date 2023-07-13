#!/bin/bash

# Define the characters to use in the combinations
chars="0123456789abcdefghijklmnopqrstuvwxyz"

# Define the base URL for downloading files
root_url="https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/"

# Define the target directory
local_dir="/home/roger/2023_self_regulatory_motifs/PDB/compressed_D"

# Create the target directory if it doesn't exist
mkdir -p $local_dir

# Iterate over each character in the first position
for i in $(echo $chars | fold -w1); do
  # Iterate over each character in the second position
  for j in $(echo $chars | fold -w1); do
    cd $local_dir
    # Define the subfolder
    subfolder=$i$j

    # Download the subfolder to the target directory
    wget -r -nd --no-parent --no-check-certificate $root_url/$subfolder/
  done
done

wait