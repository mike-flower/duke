#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_01_80c_2G_2h
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_job_$JOB_ID.out
#$ -e logs/duke_job_$JOB_ID.err
#$ -l h_rt=02:00:00
#$ -pe mpi 80
#$ -l mem=2G
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ============================================================================
# Duke Pipeline - Kathleen HPC Job Script
# ============================================================================
#
# SUBMISSION INSTRUCTIONS:
#
# 1. Configure duke_run.R with your data paths and parameters
#    IMPORTANT: Set threads = 80 (must match -pe smp 80)
#
# 2. Make this script executable:
#    chmod +x duke_kathleen.sh
#
# 3. Create logs directory (if it doesn't exist):
#    mkdir -p logs
#
# 4. Submit the job:
#    qsub duke_kathleen.sh
#
# 5. Monitor the job:
#    qstat -u $USER                          # Check job status
#    qstat -j <JOB_ID>                       # Detailed job info
#    tail -f logs/duke_job_<JOB_ID>.out      # Watch output log
#    tail -f logs/duke_job_<JOB_ID>.err      # Watch error log
#
# 6. Cancel job if needed:
#    qdel <JOB_ID>
#
# 7. After completion, check results:
#    ls -lh result_duke/
#    ls result_duke/module_data/*.RData
#
# ============================================================================
# RESOURCE PARAMETERS EXPLAINED:
# ============================================================================
#
# -N duke_kathleen_80c      Job name (appears in qstat)
# -wd /path/to/duke         Working directory (where job runs)
# -o logs/duke_job_$JOB_ID.out   Standard output log
# -e logs/duke_job_$JOB_ID.err   Standard error log
# -l h_rt=48:00:00          Walltime (48 hours max runtime)
# -pe smp 80                Parallel environment (80 CPU cores = 2 nodes)
# -l mem=2G                 Memory per CPU (160GB total = 80 × 2GB)
#
# ⚠️ CRITICAL KATHLEEN DIFFERENCES FROM MYRIAD:
#
# 1. NO tmpfs support (diskless nodes)
#    - DO NOT USE: -l tmpfs=XXG
#    - Duke uses Scratch for temp files instead
#
# 2. Cores must be multiples of 40 (node size)
#    - Valid: 40, 80, 120, 160
#    - Invalid: 12, 36, 50, 72
#
# 3. Memory limit: 4.8GB per core max
#    - 192GB per node / 40 cores = 4.8GB
#    - Requesting >4.8GB may cause long queue times
#
# 4. Exclusive nodes (no sharing)
#    - You get entire nodes (40 cores each)
#    - More predictable performance
#
# ADJUST RESOURCES FOR YOUR DATASET:
#
# Small-Medium dataset (≤500 samples):
#   -l h_rt=24:00:00 -pe smp 40 -l mem=4G
#   Expected: ~18-24 hours
#
# Large dataset (500-1000 samples):
#   -l h_rt=36:00:00 -pe smp 80 -l mem=2G
#   Expected: ~15-20 hours
#
# Very large dataset (1000-2000 samples):
#   -l h_rt=48:00:00 -pe smp 80 -l mem=3G
#   Expected: ~30-36 hours
#
# Massive dataset (2000+ samples):
#   -l h_rt=72:00:00 -pe smp 120 -l mem=2G
#   Expected: ~48-60 hours
#
# IMPORTANT: Match -pe smp value to threads parameter in duke_run.R!
# Example: -pe smp 80 requires threads = 80 in duke_run.R
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
echo "Duke Pipeline Job Started (Kathleen)"
echo "========================================"
echo "Job ID: $JOB_ID"
echo "Job Name: $JOB_NAME"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "Working Directory: $(pwd)"
echo ""
echo "Software Versions:"
echo "  R: $(R --version | head -1)"
echo "  samtools: $(samtools --version | head -1)"
echo "  minimap2: $(minimap2 --version)"
echo ""
echo "Resources:"
echo "  CPUs: $NSLOTS"
echo "  Memory: 2GB per CPU (Total: $((NSLOTS * 2))GB)"
echo "  Nodes: $((NSLOTS / 40)) (40 cores per node on Kathleen)"
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
REQUESTED_CORES=$NSLOTS
DUKE_THREADS=$(grep "^[[:space:]]*threads[[:space:]]*=" duke_run.R | sed 's/.*=[[:space:]]*\([0-9]*\).*/\1/')

if [ ! -z "$DUKE_THREADS" ] && [ "$DUKE_THREADS" != "$REQUESTED_CORES" ]; then
    echo "⚠️  WARNING: Thread mismatch detected!"
    echo "   Job script requests: $REQUESTED_CORES cores"
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
echo "========================================"
echo ""
echo "Output location:"
echo "  Check duke_run.R for dir_out parameter"
echo ""
echo "Logs:"
echo "  Job output: logs/duke_job_$JOB_ID.out"
echo "  Job errors: logs/duke_job_$JOB_ID.err"
echo ""
echo "Check results:"
echo "  ls result_duke/module_data/*.RData"
echo "  ls result_duke/*.html"
echo "========================================"

# Exit with Duke's exit code
exit $DUKE_EXIT_CODE

# ============================================================================
# TROUBLESHOOTING:
# ============================================================================
#
# Job fails immediately with "tmpfs" error:
#   ✗ Remove -l tmpfs= line (Kathleen doesn't support tmpfs)
#   ✓ Duke uses Scratch for temp files automatically
#
# Job fails with "invalid resource request":
#   ✗ Cores not multiple of 40
#   ✓ Use: -pe smp 40, 80, 120, or 160
#
# Job stays in queue (qw) for hours:
#   - Normal for Kathleen (15 min - 2 hours typical)
#   - Reduce cores if urgent (80 → 40)
#   - Check queue: qstat -j <JOB_ID>
#
# Thread mismatch warning:
#   Edit duke_run.R and set: threads = 80 (or whatever -pe smp value)
#
# Job fails but Duke runs:
#   cat logs/<TIMESTAMP>/<TIMESTAMP>_duke_run.log
#   ls -lh result_duke/module_data/
#
# Out of memory:
#   Increase: -l mem=3G or -l mem=4G
#   Enable: remove_intermediate = TRUE in duke_run.R
#
# Job exceeds walltime:
#   Increase: -l h_rt=72:00:00
#   Or use resume: Set resume = TRUE in duke_run.R and resubmit
#
# Module not found errors:
#   Check: module purge; module load r/recommended; module load samtools
#   Verify: which R; which samtools; which minimap2
#
# Resume from failure:
#   Set resume = TRUE in duke_run.R
#   Delete broken cache: rm result_duke/temp/*.RData
#   Resubmit: qsub duke_kathleen.sh
#
# Backup results (Kathleen has backed-up ACFS storage):
#   mkdir -p ~/ACFS/duke_results_$(date +%Y%m%d)
#   cp -r result_duke ~/ACFS/duke_results_$(date +%Y%m%d)/
#
# ============================================================================
# PERFORMANCE EXPECTATIONS (1000 samples):
# ============================================================================
#
# Kathleen 40 cores:  ~28 hours total
# Kathleen 80 cores:  ~15 hours total (RECOMMENDED)
# Kathleen 120 cores: ~12 hours total (diminishing returns)
#
# Module breakdown (80 cores):
#   Module 1 (Import):       ~1 hour
#   Module 2 (Alignment):    ~8 hours
#   Module 3 (Repeats):      ~5 minutes (highly parallel!)
#   Module 4 (Alleles):      ~6 hours
#   Modules 5-7 (Viz):       ~30 minutes
#
# Compare to Myriad 36 cores: ~18 hours
# Kathleen 80 cores is ~20% faster for large datasets
#
# ============================================================================
