#!/bin/bash -l
#SBATCH --job-name=duke_cli
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=48:00:00
#SBATCH --output=logs/duke_%j.out
#SBATCH --error=logs/duke_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=michael.flower@ucl.ac.uk

# ==============================================================================
# Duke Pipeline - New Kathleen (Slurm) Job Script - CLI Version
# ==============================================================================
#
# This script runs Duke using the command-line interface (duke_cli.R)
# Edit the parameters in the ./duke command below
#
# Cluster: New Kathleen (UCL)
# Scheduler: Slurm
# Cores: 80 (2 nodes × 40 cores)
# Memory: 2GB per core = 160GB total
# Runtime: 48 hours
#
# CRITICAL NEW KATHLEEN (SLURM) DIFFERENCES:
# 1. Uses Slurm scheduler (NOT SGE)
#    - Submit: sbatch (not qsub)
#    - Check: squeue (not qstat)
#    - Cancel: scancel (not qdel)
#    - NO chmod needed: Slurm scripts don't need to be executable
#
# 2. Minimum 80 cores (2 nodes × 40 cores each)
#    - Valid: 80, 120, 160 (multiples of 40)
#    - Invalid: 12, 36, 40, 72
#
# 3. NO tmpfs support (diskless nodes)
#    - Duke uses Scratch for temp files automatically
#
# 4. Exclusive nodes (no sharing)
#    - More predictable performance
#
# Expected performance (80 cores):
#   ~1000 samples: 15-18 hours
#   ~2000 samples: 30-36 hours
#
# ==============================================================================

# ==============================================================================
# SUBMISSION INSTRUCTIONS:
# ==============================================================================
#
# 1. Edit parameters in the ./duke command below (see "EDIT PARAMETERS" section)
#
# 2. Create logs directory (if it doesn't exist):
#    mkdir -p logs
#
# 3. Submit the job:
#    sbatch duke_kathleen_slurm_cli.sh
#
#    Note: NO chmod needed for Slurm scripts (unlike SGE)
#
# 4. Monitor the job:
#    squeue -u $USER                         # Check job status
#    scontrol show job <JOB_ID>              # Detailed job info
#    tail -f logs/duke_<JOB_ID>.out          # Watch output log
#
# 5. Cancel job if needed:
#    scancel <JOB_ID>
#
# ==============================================================================
# RESOURCE PARAMETERS EXPLAINED:
# ==============================================================================
#
# --job-name=duke_cli       Job name (appears in squeue)
# --nodes=2                 Number of compute nodes (2 × 40 cores = 80 total)
# --ntasks=80               Total cores (must match --threads parameter!)
# --cpus-per-task=1         CPUs per task (always 1 for Duke)
# --mem-per-cpu=2G          Memory per CPU (80 × 2G = 160GB total)
# --time=48:00:00           Walltime in HH:MM:SS format
# --output=logs/duke_%j.out Standard output (%j = job ID)
# --error=logs/duke_%j.err  Standard error (%j = job ID)
#
# ADJUST RESOURCES FOR YOUR DATASET:
#
# Small-Medium (≤500 samples):
#   --time=24:00:00 --ntasks=80 --mem-per-cpu=4G
#
# Large (500-1000 samples):
#   --time=36:00:00 --ntasks=80 --mem-per-cpu=2G
#
# Very large (1000-2000 samples):
#   --time=48:00:00 --ntasks=80 --mem-per-cpu=3G
#
# Massive (2000+ samples):
#   --time=72:00:00 --ntasks=120 --mem-per-cpu=2G --nodes=3
#
# IMPORTANT: Match --ntasks to --threads in duke command!
# Example: --ntasks=80 requires --threads 80
#
# ==============================================================================

# ==============================================================================
# Load Modules
# ==============================================================================
module purge
module load r/recommended
module load samtools/1.11/gnu-4.9.2

# Set R library path
export R_LIBS_USER=~/R/library

# Add minimap2 to PATH
export PATH=$HOME/Scratch/bin/minimap2:$PATH

# ==============================================================================
# Print Environment Info
# ==============================================================================
echo "========================================"
echo "Duke Pipeline - New Kathleen (Slurm CLI)"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks (cores): $SLURM_NTASKS"
echo "Memory: $((SLURM_NTASKS * 2))GB total"
echo "Walltime: 48 hours"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo ""
echo "Software Versions:"
echo "  R: $(R --version | head -1)"
echo "  samtools: $(samtools --version | head -1)"
echo "  minimap2: $(minimap2 --version)"
echo "========================================"
echo ""

# ==============================================================================
# Pre-flight Checks
# ==============================================================================
echo "Running pre-flight checks..."

# Check duke_cli.R exists
if [ ! -f "duke" ]; then
    echo "ERROR: duke CLI wrapper not found in $(pwd)"
    exit 1
fi
echo "✓ duke CLI found"

# Check lib directory exists
if [ ! -d "lib" ]; then
    echo "ERROR: lib/ directory not found in $(pwd)"
    exit 1
fi
echo "✓ lib/ directory found"

# Check minimap2 is available
if ! command -v minimap2 &> /dev/null; then
    echo "ERROR: minimap2 not found in PATH"
    exit 1
fi
echo "✓ minimap2 found: $(which minimap2)"

echo ""
echo "Pre-flight checks complete. Starting Duke Pipeline..."
echo ""

# ==============================================================================
# EDIT PARAMETERS BELOW
# ==============================================================================

# Change to Duke directory
cd ~/Scratch/bin/duke

# Run Duke with CLI
# IMPORTANT: Keep --threads 80 to match --ntasks=80 above
# Edit these paths for your analysis:

./duke \
  --dir_data /home/skgtmdf/Scratch/data/2025.12.17_pb_test/data \
  --dir_out /home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke \
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta \
  --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads 80 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE

# Capture exit code
DUKE_EXIT_CODE=$?

# ==============================================================================
# COMMON CUSTOMISATIONS
# ==============================================================================

# Without adapter trimming:
# ./duke \
#   --dir_data ~/Scratch/data/my_experiment \
#   --dir_out ~/Scratch/results/my_results \
#   --path_ref ~/Scratch/refs/HTTset20/HTTset20.fasta \
#   --trim FALSE \
#   --threads 80

# Force re-run all modules:
# ./duke --resume FALSE ...

# Run specific modules only (e.g., re-plot after parameter change):
# ./duke --run_modules 5,6,7 --resume TRUE ...

# Downsample for testing:
# ./duke --downsample 1000 ...

# Custom repeat parameters:
# ./duke --rpt_pattern CTG --rpt_max_mismatch 1 ...

# ==============================================================================
# Print Completion Message
# ==============================================================================
echo ""
echo "========================================"
if [ $DUKE_EXIT_CODE -eq 0 ]; then
    echo "✓ Duke Pipeline Job Completed Successfully"
else
    echo "✗ Duke Pipeline Job Failed (exit code: $DUKE_EXIT_CODE)"
fi
echo "Date: $(date)"
echo "Duration: $(date -ud "@$SECONDS" +%T)"
echo "========================================"
echo ""
echo "Logs:"
echo "  Job output: logs/duke_${SLURM_JOB_ID}.out"
echo "  Job errors: logs/duke_${SLURM_JOB_ID}.err"
echo ""
echo "Check results:"
echo "  ls result_duke/module_data/*.RData"
echo "  ls result_duke/*.html"
echo "========================================"

# Exit with Duke's exit code
exit $DUKE_EXIT_CODE

# ==============================================================================
# SLURM COMMANDS REFERENCE:
# ==============================================================================
#
# Submit job:
#   sbatch duke_kathleen_slurm_cli.sh
#
# Check status:
#   squeue -u $USER
#   squeue -j <JOB_ID>
#
# Job details:
#   scontrol show job <JOB_ID>
#
# Cancel job:
#   scancel <JOB_ID>
#
# Watch queue (updates every 30 seconds):
#   watch -n 30 'squeue -u $USER'
#
# After completion:
#   sacct -j <JOB_ID>
#   sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,MaxRSS
#
# ==============================================================================
# TROUBLESHOOTING:
# ==============================================================================
#
# Job stays in pending (PD) state:
#   - Normal: Queue times vary (15 min - 2 hours)
#   - Check reason: squeue -j <JOB_ID>
#
# Job fails immediately:
#   - Check: cat logs/duke_<JOB_ID>.err
#   - Common causes:
#     * Module not found (check module load commands)
#     * R packages missing
#     * minimap2 not in PATH
#
# Thread mismatch:
#   Make sure --threads in duke command matches --ntasks in header
#   Example: --ntasks=80 requires --threads 80
#
# Out of memory:
#   Increase: --mem-per-cpu=3G or --mem-per-cpu=4G
#   Enable: --remove_intermediate TRUE
#
# Job exceeds walltime:
#   Increase: --time=72:00:00
#   Or use resume: --resume TRUE and resubmit
#
# Resume from failure:
#   ./duke command with --resume TRUE
#   Delete broken cache: rm result_duke/module_data/<failed_module>.RData
#   Resubmit: sbatch duke_kathleen_slurm_cli.sh
#
# ==============================================================================
# SLURM VS SGE COMMAND COMPARISON:
# ==============================================================================
#
# | Action         | SGE (Old Kathleen) | Slurm (New Kathleen)        |
# |----------------|--------------------|-----------------------------|
# | Make executable| chmod +x script.sh | NOT NEEDED                  |
# | Submit job     | qsub script.sh     | sbatch script.sh            |
# | Check jobs     | qstat -u $USER     | squeue -u $USER             |
# | Job details    | qstat -j JOB_ID    | scontrol show job JOB_ID    |
# | Cancel job     | qdel JOB_ID        | scancel JOB_ID              |
# | After complete | qacct -j JOB_ID    | sacct -j JOB_ID             |
# | Cluster info   | qhost -q           | sinfo                       |
# | Job ID var     | $JOB_ID            | $SLURM_JOB_ID               |
#
# ==============================================================================
# PERFORMANCE EXPECTATIONS:
# ==============================================================================
#
# Test run (3-10 samples, 80 cores):
#   3 samples:  ~20 minutes
#   10 samples: ~55 minutes
#
# Production (1000 samples, 80 cores):
#   Module 1 (Import):       ~1 hour
#   Module 2 (Alignment):    ~8 hours
#   Module 3 (Repeats):      ~5 minutes
#   Module 4 (Alleles):      ~6 hours
#   Modules 5-7 (Viz):       ~30 minutes
#   Total:                   ~15 hours
#
# Compare to Old Kathleen (SGE):
#   Same performance (same hardware)
#   New Kathleen may have shorter queue times
#
# Compare to Myriad:
#   Myriad 36 cores: ~52 hours (1000 samples)
#   Kathleen 80 cores: ~15 hours (3.5x faster)
#
# ==============================================================================
