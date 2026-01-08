#!/bin/bash -l
#SBATCH --job-name=duke_test
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00
#SBATCH --output=logs/duke_%j.out
#SBATCH --error=logs/duke_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=michael.flower@ucl.ac.uk

# ==============================================================================
# Duke Pipeline - New Kathleen (Slurm) - Test Run
# ==============================================================================
# Quick test run with 80 cores, 2G memory, 2 hours
# For small test datasets (3-10 samples)
# ==============================================================================

module purge
module load r/recommended
module load samtools/1.11/gnu-4.9.2

export R_LIBS_USER=~/R/library
export PATH=$HOME/Scratch/bin/minimap2:$PATH

# ==============================================================================
# Print Environment Info
# ==============================================================================
echo "========================================"
echo "Duke Pipeline - New Kathleen (Slurm)"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Memory per CPU: 2G (160G total)"
echo "Walltime: 2 hours"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "Working Directory: $(pwd)"
echo ""
echo "Software Versions:"
echo "  R: $(R --version | head -1)"
echo "  samtools: $(samtools --version | head -1)"
echo "  minimap2: $(minimap2 --version)"
echo ""
echo "========================================"
echo ""

# ==============================================================================
# Pre-flight Checks
# ==============================================================================
echo "Running pre-flight checks..."

if [ ! -f "duke_run.R" ]; then
    echo "ERROR: duke_run.R not found in $(pwd)"
    exit 1
fi
echo "✓ duke_run.R found"

if [ ! -d "lib" ]; then
    echo "ERROR: lib/ directory not found in $(pwd)"
    exit 1
fi
echo "✓ lib/ directory found"

if ! command -v minimap2 &> /dev/null; then
    echo "ERROR: minimap2 not found in PATH"
    exit 1
fi
echo "✓ minimap2 found: $(which minimap2)"

# Check threads match
REQUESTED_TASKS=$SLURM_NTASKS
DUKE_THREADS=$(grep "^[[:space:]]*threads[[:space:]]*=" duke_run.R | sed 's/.*=[[:space:]]*\([0-9]*\).*/\1/')

if [ ! -z "$DUKE_THREADS" ] && [ "$DUKE_THREADS" != "$REQUESTED_TASKS" ]; then
    echo "⚠️  WARNING: Thread mismatch!"
    echo "   Job requests: $REQUESTED_TASKS cores"
    echo "   duke_run.R:   $DUKE_THREADS cores"
    echo "   These should match for optimal performance."
fi

echo ""
echo "Pre-flight checks complete. Starting Duke Pipeline..."
echo ""

# ==============================================================================
# Run Duke Pipeline
# ==============================================================================
cd ~/Scratch/bin/duke
Rscript duke_run.R

DUKE_EXIT_CODE=$?

# ==============================================================================
# Completion Message
# ==============================================================================
echo ""
echo "========================================"
if [ $DUKE_EXIT_CODE -eq 0 ]; then
    echo "✓ Duke Pipeline Completed Successfully"
else
    echo "✗ Duke Pipeline Failed (exit code: $DUKE_EXIT_CODE)"
fi
echo "Date: $(date)"
echo "Duration: $(date -ud "@$SECONDS" +%T)"
echo "========================================"
echo ""
echo "Results location: Check duke_run.R for dir_out"
echo "Logs: logs/duke_${SLURM_JOB_ID}.out"
echo ""
echo "Check results:"
echo "  ls result_duke/module_data/*.RData"
echo "  ls result_duke/*.html"
echo "========================================"

exit $DUKE_EXIT_CODE

# ==============================================================================
# SLURM COMMANDS REFERENCE
# ==============================================================================
#
# Submit job:
#   sbatch duke_kathleen_slurm.sh
#
# Check status:
#   squeue -u $USER
#
# Job details:
#   scontrol show job JOB_ID
#
# Cancel job:
#   scancel JOB_ID
#
# Watch queue:
#   watch -n 30 'squeue -u $USER'
#
# After completion:
#   sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,MaxRSS
#
# ==============================================================================
# EXPECTED PERFORMANCE (test run, 3-10 samples, 80 cores)
# ==============================================================================
#
# 3 samples:
#   Module 1:  ~2 min
#   Module 2:  ~10 min
#   Module 3:  ~5 sec
#   Module 4:  ~5 min
#   Modules 5-7: ~2 min
#   Total:     ~20 min
#
# 10 samples:
#   Module 1:  ~5 min
#   Module 2:  ~30 min
#   Module 3:  ~10 sec
#   Module 4:  ~15 min
#   Modules 5-7: ~5 min
#   Total:     ~55 min
#
# 2 hour walltime provides plenty of buffer for test runs
#
# ==============================================================================
