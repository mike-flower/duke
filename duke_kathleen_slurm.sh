#!/bin/bash -l
#SBATCH --job-name=duke_01_2n_80c_2G_2h
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00
#SBATCH --output=logs/duke_%j.out
#SBATCH --error=logs/duke_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=michael.flower@ucl.ac.uk

# ============================================================================
# Duke Pipeline - New Kathleen (Slurm) Job Script
# ============================================================================
#
# SUBMISSION INSTRUCTIONS:
#
# 1. Configure duke_run.R with your data paths and parameters
#    IMPORTANT: Set threads = 80 (must match --ntasks=80)
#
# 2. Create logs directory (if it doesn't exist):
#    mkdir -p logs
#
# 3. Submit the job:
#    sbatch duke_kathleen_slurm.sh
#
#    Note: Slurm scripts do NOT need chmod +x (unlike SGE scripts)
#
# 4. Monitor the job:
#    squeue -u $USER                         # Check job status
#    scontrol show job <JOB_ID>              # Detailed job info
#    tail -f logs/duke_<JOB_ID>.out          # Watch output log
#    tail -f logs/duke_<JOB_ID>.err          # Watch error log
#
# 5. Cancel job if needed:
#    scancel <JOB_ID>
#
# 6. After completion, check results:
#    ls -lh result_duke/
#    ls result_duke/module_data/*.RData
#
# ============================================================================
# RESOURCE PARAMETERS EXPLAINED:
# ============================================================================
#
# --job-name=duke_test      Job name (appears in squeue)
# --nodes=2                 Number of compute nodes (2 × 40 cores = 80 total)
# --ntasks=80               Total number of tasks/cores (must match threads in duke_run.R!)
# --cpus-per-task=1         CPUs per task (always 1 for Duke)
# --mem-per-cpu=2G          Memory per CPU (80 × 2G = 160GB total)
# --time=2:00:00            Walltime in HH:MM:SS format (2 hours)
# --output=logs/duke_%j.out Standard output (%j = job ID)
# --error=logs/duke_%j.err  Standard error (%j = job ID)
# --mail-type=...           When to send email (BEGIN,END,FAIL)
# --mail-user=...           Email address for notifications
#
# ⚠️ CRITICAL NEW KATHLEEN (SLURM) DIFFERENCES:
#
# 1. Uses Slurm scheduler (NOT SGE)
#    - Submit: sbatch (not qsub)
#    - Check: squeue (not qstat)
#    - Cancel: scancel (not qdel)
#    - NO chmod needed: Slurm scripts don't need to be executable
#
# 2. Script permissions:
#    - SGE: Requires chmod +x script.sh before qsub
#    - Slurm: NO chmod needed, sbatch works on any .sh file
#
# 3. Minimum 80 cores (2 nodes × 40 cores each)
# 3. Minimum 80 cores (2 nodes × 40 cores each)
#    - Valid: 80, 120, 160 (multiples of 40)
#    - Invalid: 12, 36, 40, 72
#
# 4. NO tmpfs support (diskless nodes)
#    - DO NOT USE: tmpfs or local storage
#    - Duke uses Scratch for temp files automatically
#
# 5. Exclusive nodes (no sharing)
#    - You get entire nodes
#    - More predictable performance
#
# ADJUST RESOURCES FOR YOUR DATASET:
#
# Small-Medium dataset (≤500 samples):
#   --time=24:00:00 --ntasks=80 --mem-per-cpu=4G
#   Expected: ~18-24 hours
#
# Large dataset (500-1000 samples):
#   --time=36:00:00 --ntasks=80 --mem-per-cpu=2G
#   Expected: ~15-20 hours
#
# Very large dataset (1000-2000 samples):
#   --time=48:00:00 --ntasks=80 --mem-per-cpu=3G
#   Expected: ~30-36 hours
#
# Massive dataset (2000+ samples):
#   --time=72:00:00 --ntasks=120 --mem-per-cpu=2G --nodes=3
#   Expected: ~48-60 hours
#
# IMPORTANT: Match --ntasks value to threads parameter in duke_run.R!
# Example: --ntasks=80 requires threads = 80 in duke_run.R
#
# ============================================================================

# ============================================================================
# Load Modules (ALL software dependencies)
# ============================================================================
module purge
module load r/recommended                    # R 4.2.0
module load samtools/1.11/gnu-4.9.2         # Samtools 1.11

# Set R library path
export R_LIBS_USER=~/R/library

# Add minimap2 (compiled in ~/Scratch/bin/minimap2)
export PATH=$HOME/Scratch/bin/minimap2:$PATH

# ============================================================================
# Print Environment Info
# ============================================================================
echo "========================================"
echo "Duke Pipeline - New Kathleen (Slurm)"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks (cores): $SLURM_NTASKS"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Memory per CPU: 2G (Total: $((SLURM_NTASKS * 2))GB)"
echo "Walltime limit: 2 hours"
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

# ============================================================================
# Pre-flight Checks
# ============================================================================
echo "Running pre-flight checks..."

# Check duke_run.R exists
if [ ! -f "duke_run.R" ]; then
    echo "ERROR: duke_run.R not found in $(pwd)"
    exit 1
fi
echo "✓ duke_run.R found"

# Check lib directory exists
if [ ! -d "lib" ]; then
    echo "ERROR: lib/ directory not found in $(pwd)"
    exit 1
fi
echo "✓ lib/ directory found"

# Check minimap2 is available
if ! command -v minimap2 &> /dev/null; then
    echo "ERROR: minimap2 not found in PATH"
    echo "Expected at: $HOME/Scratch/bin/minimap2/minimap2"
    exit 1
fi
echo "✓ minimap2 found: $(which minimap2)"

# Check threads parameter matches job request
REQUESTED_TASKS=$SLURM_NTASKS
DUKE_THREADS=$(grep "^[[:space:]]*threads[[:space:]]*=" duke_run.R | sed 's/.*=[[:space:]]*\([0-9]*\).*/\1/')

if [ ! -z "$DUKE_THREADS" ] && [ "$DUKE_THREADS" != "$REQUESTED_TASKS" ]; then
    echo "⚠️  WARNING: Thread mismatch detected!"
    echo "   Job script requests: $REQUESTED_TASKS cores"
    echo "   duke_run.R threads:  $DUKE_THREADS"
    echo "   For optimal performance, these should match."
    echo "   Continuing anyway..."
fi

echo ""
echo "Pre-flight checks complete. Starting Duke Pipeline..."
echo ""

# ============================================================================
# Run Duke Pipeline
# ============================================================================
cd ~/Scratch/bin/duke
Rscript duke_run.R

# Capture exit code
DUKE_EXIT_CODE=$?

# ============================================================================
# Print Completion Message
# ============================================================================
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
echo "Output location:"
echo "  Check duke_run.R for dir_out parameter"
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

# ============================================================================
# SLURM COMMANDS REFERENCE:
# ============================================================================
#
# Submit job:
#   sbatch duke_kathleen_slurm.sh
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
# Check cluster info:
#   sinfo
#   sinfo -N -l
#
# ============================================================================
# TROUBLESHOOTING:
# ============================================================================
#
# Job stays in pending (PD) state:
#   - Normal: Queue times vary (15 min - 2 hours typical)
#   - Check reason: squeue -j <JOB_ID>
#   - Reduce cores if urgent (80 → contact RC support for options)
#
# Job fails immediately:
#   - Check: cat logs/duke_<JOB_ID>.err
#   - Common causes:
#     * Module not found (check module purge; module load commands)
#     * R packages missing (re-run R package installation)
#     * minimap2 not in PATH (check export PATH line)
#
# Thread mismatch warning:
#   Edit duke_run.R and set: threads = 80 (or whatever --ntasks value)
#
# Job fails but Duke runs:
#   cat logs/<TIMESTAMP>/<TIMESTAMP>_duke_run.log
#   ls -lh result_duke/module_data/
#
# Out of memory:
#   Increase: --mem-per-cpu=3G or --mem-per-cpu=4G
#   Enable: remove_intermediate = TRUE in duke_run.R
#
# Job exceeds walltime:
#   Increase: --time=72:00:00
#   Or use resume: Set resume = TRUE in duke_run.R and resubmit
#
# Module not found errors:
#   Check: module purge; module load r/recommended; module load samtools
#   Verify: which R; which samtools; which minimap2
#
# Resume from failure:
#   Set resume = TRUE in duke_run.R
#   Delete broken cache: rm result_duke/temp/*.RData
#   Resubmit: sbatch duke_kathleen_slurm.sh
#
# Backup results (Kathleen has backed-up ACFS storage):
#   mkdir -p ~/ACFS/duke_results_$(date +%Y%m%d)
#   cp -r result_duke ~/ACFS/duke_results_$(date +%Y%m%d)/
#
# Cannot submit to Slurm on New Kathleen:
#   - New Kathleen might not have Slurm yet (check: sinfo)
#   - Use Old Kathleen with SGE instead: ssh kathleen.rc.ucl.ac.uk
#   - Use qsub with duke_kathleen.sh (SGE version)
#
# ============================================================================
# PERFORMANCE EXPECTATIONS:
# ============================================================================
#
# Test run (3-10 samples, 80 cores):
#   3 samples:  ~20 minutes total
#   10 samples: ~55 minutes total
#
# Production (1000 samples, 80 cores):
#   Module 1 (Import):       ~1 hour
#   Module 2 (Alignment):    ~8 hours
#   Module 3 (Repeats):      ~5 minutes (highly parallel!)
#   Module 4 (Alleles):      ~6 hours
#   Modules 5-7 (Viz):       ~30 minutes
#   Total:                   ~15 hours
#
# Compare different core counts (1000 samples):
#   80 cores:  ~15 hours total (RECOMMENDED)
#   120 cores: ~12 hours total (diminishing returns)
#
# Compare to Old Kathleen (SGE):
#   Same performance expected (both are Kathleen hardware)
#   New Kathleen may have shorter queue times (less crowded)
#
# Compare to Myriad:
#   Myriad 36 cores: ~18 hours (1000 samples)
#   Kathleen 80 cores: ~15 hours (20% faster)
#
# ============================================================================
# SLURM VS SGE COMMAND COMPARISON:
# ============================================================================
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
# ============================================================================
