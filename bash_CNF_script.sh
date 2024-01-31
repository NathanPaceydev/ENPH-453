#!/bin/bash

# Check if there are any .CNF files in the current directory
shopt -s nullglob
cnf_files=(*.CNF)

if [ ${#cnf_files[@]} -eq 0 ]; then
    echo "No .CNF files found in the current directory."
    exit 1
fi

# Loop through every .CNF file in the current directory
for file in "${cnf_files[@]}"; do
    echo "Processing file: $file"
    python3 particleCNFconverter.py "$file"
    
    #echo "Plotting file: $file"
    #python3 process_cnf.py "$file.txt"

done
