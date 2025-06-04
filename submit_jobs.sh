#!/bin/bash

# Configuration
SEEDS=(1 2 3 4 5 6 7 8 9 10)  # Add or modify seed values as needed
BASE_PBS_FILE="run_sim.pbs"
TEMP_PBS_PREFIX="run_sim_seed"

# Check if base PBS file exists
if [ ! -f "$BASE_PBS_FILE" ]; then
    echo "Error: Base PBS file '$BASE_PBS_FILE' not found!"
    exit 1
fi

echo "Submitting jobs with seeds: ${SEEDS[@]}"

# Loop through each seed value
for seed in "${SEEDS[@]}"; do
    # Create temporary PBS file name
    temp_pbs_file="${TEMP_PBS_PREFIX}_${seed}.pbs"
    
    # Copy base PBS file and modify it
    cp "$BASE_PBS_FILE" "$temp_pbs_file"
    
    # Replace the main command line to include seed and update name
    sed -i.bak "s|--name smc|--seed $seed --name smc_seed_$seed|g" "$temp_pbs_file"
    
    # Remove backup file
    rm "${temp_pbs_file}.bak"
    
    echo "Submitting job with seed $seed..."
    
    # Submit the job
    # qsub "$temp_pbs_file"
    
    # Check if submission was successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully submitted job with seed $seed"
    else
        echo "✗ Failed to submit job with seed $seed"
    fi
    
    # Clean up temporary PBS file
    # rm "$temp_pbs_file"
    
    # Optional: Add a small delay between submissions
    sleep 1
done

echo "All jobs submitted!" 