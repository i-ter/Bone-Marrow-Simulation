#!/bin/bash

# Default values
N_JOBS=0
NAME=""
TEMP_PBS_PREFIX="run"
DRY_RUN=false

# Define the swap_motility values to loop over
SWAP_MOTILITY_VALUES=( 0.05 0.10 0.20 0.30 0.50 1.00)

usage() {
    echo "Usage: $0 --name <job_name> -n_jobs <number> [--dry-run]"
    echo "  --name     Name for the job"
    echo "  -n_jobs    Number of jobs to submit per swap_motility value (default: 10)"
    echo "  --dry-run  Show what jobs would be submitted without actually submitting them"
    echo ""
    echo "This script will submit jobs for each swap_motility value: ${SWAP_MOTILITY_VALUES[*]}"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --name)
            NAME="$2"
            shift 2
            ;;
        -n_jobs)
            N_JOBS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
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

if ! [[ "$N_JOBS" =~ ^[0-9]+$ ]] || [ "$N_JOBS" -le 0 ]; then
    echo "Error: -n_jobs must be a positive integer"
    usage
fi

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN MODE - No jobs will actually be submitted"
fi
echo "Submitting $N_JOBS jobs per swap_motility value with base name: $NAME"
echo "Swap motility values: ${SWAP_MOTILITY_VALUES[*]}"

total_jobs=$((N_JOBS * ${#SWAP_MOTILITY_VALUES[@]}))
echo "Total jobs to submit: $total_jobs"

# Loop through each swap_motility value
for swap_motility in "${SWAP_MOTILITY_VALUES[@]}"; do
    echo ""
    echo "Submitting jobs for swap_motility = $swap_motility"
    
    # Convert swap_motility to a filename-safe string (remove .)
    swap_motility_str=$(echo "$swap_motility" | sed 's/\.//g')
    
    # Loop through each job for this swap_motility value
    for ((i=1; i<=N_JOBS; i++)); do
        temp_pbs_file="${TEMP_PBS_PREFIX}_swapm_${swap_motility_str}_seed_${i}.pbs"
        job_name="${NAME}_swapm_${swap_motility_str}_seed_$i"
        
        # Create PBS file from scratch
        cat > "$temp_pbs_file" << EOF
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=20gb:ompthreads=12:cpu_type=icelake

ml GCC

\$PBS_O_WORKDIR/main --width 500 --height 500 --steps 200000 --cold_start true --cells 3 --swap_motility $swap_motility -dwf 500 --seed $i --name $job_name
EOF
        echo "  Submitting job $i (seed $i) for swap_motility $swap_motility..."
        if [ "$DRY_RUN" = true ]; then
            echo "    DRY RUN: Would execute: qsub $temp_pbs_file"
            echo "    PBS file content:"
            cat "$temp_pbs_file" | sed 's/^/      /'
        else
            qsub "$temp_pbs_file"
        fi
        sleep 1
    done
done

echo ""
echo "All $total_jobs jobs submitted!"
echo "Job naming convention: ${NAME}_sm[VALUE]_seed_[SEED]"
echo "  where [VALUE] is swap_motility with dots replaced by underscores"
echo "  and [SEED] is the random seed number"