#!/bin/bash

# Default values
N_JOBS=0
NAME=""
SMC=""
BASE_PBS_FILE="run_sim.pbs"
TEMP_PBS_PREFIX="run_sim_seed"

# Function to display usage
usage() {
    echo "Usage: $0 --name <job_name> --smc <true/false> -n_jobs <number>"
    echo "  --name     Name for the job"
    echo "  --smc      stop motility flag (true or false)"
    echo "  -n_jobs    Number of jobs to submit (default: 10)"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --name)
            NAME="$2"
            shift 2
            ;;
        --smc)
            SMC="$2"
            shift 2
            ;;
        -n_jobs)
            N_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$NAME" ]; then
    echo "Error: --name is required"
    usage
fi

if [ -z "$SMC" ]; then
    echo "Error: --smc is required"
    usage
fi

if [[ "$SMC" != "true" && "$SMC" != "false" ]]; then
    echo "Error: --smc must be 'true' or 'false'"
    usage
fi

if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]] || [ "$N_JOBS" -le 0 ]; then
    echo "Error: -n_jobs must be a positive integer"
    usage
fi

# Check if base PBS file exists
if [ ! -f "$BASE_PBS_FILE" ]; then
    echo "Error: Base PBS file '$BASE_PBS_FILE' not found!"
    exit 1
fi

echo "Submitting $N_JOBS jobs with name: $NAME, smc: $SMC"

# Loop through each job
for ((i=1; i<=N_JOBS; i++)); do
    # Create temporary PBS file name
    temp_pbs_file="${TEMP_PBS_PREFIX}_${i}.pbs"
    
    # Create PBS file from scratch with correct arguments
    cat > "$temp_pbs_file" << EOF
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=20gb:ompthreads=12:cpu_type=icelake

ml GCC

\$PBS_O_WORKDIR/main --width 500 --height 500 --steps 300000 --cold_start true --cells 3 -smc $SMC -dwf 500 --seed $i --name ${NAME}_seed_$i
EOF
    
    echo "Submitting job $i with seed $i..."
    
    # Submit the job
    # qsub "$temp_pbs_file"
    
    # Clean up temporary PBS file
    # rm "$temp_pbs_file"
    
    # Optional: Add a small delay between submissions
    sleep 1
done

echo "All $N_JOBS jobs submitted!"