#!/bin/bash

# Default values
N_JOBS=0
NAME=""
TEMP_PBS_PREFIX="run"
DRY_RUN=false
SEED_BASELINE=100

# Define the swap_motility values to loop over
SWAP_MOTILITY_VALUES=(0.00 0.05 0.10 0.20 0.30 0.50 1.00)

usage() {
    echo "Usage: $0 --name <job_name> -n_jobs <number> [--dry-run]"
    echo "  --name     Name for the job array"
    echo "  -n_jobs    Number of jobs to submit per swap_motility value"
    echo "  --dry-run  Show what job would be submitted without actually submitting it"
    echo ""
    echo "This script will submit an array job for each swap_motility value: ${SWAP_MOTILITY_VALUES[*]}"
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

total_jobs=$((N_JOBS * ${#SWAP_MOTILITY_VALUES[@]}))

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN MODE - No job will actually be submitted"
fi
echo "Preparing to submit an array job of size $total_jobs for base name: $NAME"
echo "This comprises $N_JOBS jobs for each of the following swap motility values: ${SWAP_MOTILITY_VALUES[*]}"


temp_pbs_file="array.pbs"

# Create PBS file for the array job
cat > "$temp_pbs_file" << EOF
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=20gb:ompthreads=12:cpu_type=icelake
#PBS -J 1-${total_jobs}
#PBS -o ${NAME}_\$PBS_ARRAY_INDEX.out
#PBS -e ${NAME}_\$PBS_ARRAY_INDEX.err

# ---BEGIN VARIABLES---
# These variables are set by submit_jobs.sh
SWAP_MOTILITY_VALUES=(${SWAP_MOTILITY_VALUES[*]})
N_JOBS=${N_JOBS}
SEED_BASELINE=${SEED_BASELINE}
NAME="${NAME}"
# ---END VARIABLES---

# Calculate parameters for this sub-job based on PBS_ARRAY_INDEX
task_id=\$((PBS_ARRAY_INDEX - 1))
swap_idx=\$((task_id / N_JOBS))
job_idx=\$((task_id % N_JOBS))

swap_motility=\${SWAP_MOTILITY_VALUES[\$swap_idx]}
seed=\$((SEED_BASELINE + job_idx))

# Create a filename-safe string for swap_motility
swap_motility_str=\$(echo "\$swap_motility" | sed 's/\.//g')

# Construct job name for this sub-job
job_name="\${NAME}_\${PBS_ARRAY_INDEX}_swapm_\${swap_motility_str}_seed_\${seed}"

# Load required modules
ml GCC

# Execute the simulation
\$PBS_O_WORKDIR/main --width 500 --height 500 --steps 200000 --cold_start true --cells 3 --swap_motility \$swap_motility -dwf 500 --seed \$seed --name \$job_name
EOF

echo ""
echo "Generated PBS array job script: $temp_pbs_file"

if [ "$DRY_RUN" = true ]; then
    echo "    DRY RUN: Would execute: qsub $temp_pbs_file"
    echo "--- Start of ${temp_pbs_file} ---"
    cat "$temp_pbs_file"
    echo "--- End of ${temp_pbs_file} ---"
else
    echo "Submitting array job..."
    qsub "$temp_pbs_file"
fi

echo ""
echo "Array job submission script finished."