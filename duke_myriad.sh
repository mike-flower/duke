#!/bin/bash -l
#$ -N duke_pipeline
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_job_$JOB_ID.out
#$ -e logs/duke_job_$JOB_ID.err
#$ -l h_rt=24:00:00
#$ -pe smp 4
#$ -l mem=8G
#$ -l tmpfs=50G

# ============================================================================
# Duke Pipeline - Myriad HPC Job Script
# ============================================================================
#
# SUBMISSION INSTRUCTIONS:
#
# 1. Configure duke_run_myriad.R with your data paths and parameters
#
# 2. Make this script executable:
#    chmod +x duke_myriad.sh
#
# 3. Create logs directory (if it doesn't exist):
#    mkdir -p logs
#
# 4. Submit the job:
#    qsub duke_myriad.sh
#
# 5. Monitor the job:
#    qstat                                    # Check job status
#    qstat -j <JOB_ID>                        # Detailed job info
#    tail -f logs/duke_job_<JOB_ID>.out      # Watch output log
#    tail -f logs/duke_job_<JOB_ID>.err      # Watch error log
#
# 6. Cancel job if needed:
#    qdel <JOB_ID>
#
# 7. After completion, check results:
#    ls -lh ~/Scratch/data/YOUR_PROJECT/result_duke/
#    tree ~/Scratch/data/YOUR_PROJECT/result_duke/
#
# ============================================================================
# RESOURCE PARAMETERS EXPLAINED:
# ============================================================================
#
# -N duke_pipeline          Job name (appears in qstat)
# -wd /path/to/duke         Working directory (where job runs)
# -o logs/duke_job_$JOB_ID.out   Standard output log
# -e logs/duke_job_$JOB_ID.err   Standard error log
# -l h_rt=24:00:00          Walltime (24 hours max runtime)
# -pe smp 4                 Parallel environment (4 CPU cores)
# -l mem=8G                 Memory per CPU (32GB total = 4 × 8GB)
# -l tmpfs=50G              Temporary disk space for BAM files
#
# ADJUST RESOURCES FOR YOUR DATASET:
#
# Small dataset (1-5 samples, <50K reads/sample):
#   -l h_rt=12:00:00 -pe smp 2 -l mem=4G -l tmpfs=20G
#
# Medium dataset (5-20 samples, 50-200K reads/sample):
#   -l h_rt=24:00:00 -pe smp 4 -l mem=8G -l tmpfs=50G
#
# Large dataset (20-50 samples, 200K-1M reads/sample):
#   -l h_rt=48:00:00 -pe smp 8 -l mem=16G -l tmpfs=100G
#
# Very large dataset (50+ samples, >1M reads/sample):
#   -l h_rt=72:00:00 -pe smp 16 -l mem=32G -l tmpfs=200G
#
# IMPORTANT: Match -pe smp value to threads parameter in duke_run_myriad.R!
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
echo "Duke Pipeline Job Started"
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
echo "  Memory: 8GB per CPU (Total: $((NSLOTS * 8))GB)"
echo "  Temp space: 50GB"
echo "========================================"
echo ""

# ============================================================================
# Run Duke Pipeline
# ============================================================================
cd ~/Scratch/bin/duke
Rscript duke_run_myriad.R

# ============================================================================
# Print Completion Message
# ============================================================================
echo ""
echo "========================================"
echo "Duke Pipeline Job Completed"
echo "Date: $(date)"
echo "========================================"
echo ""
echo "Output location:"
echo "  Check duke_run_myriad.R for dir_out parameter"
echo ""
echo "Logs:"
echo "  Job output: logs/duke_job_$JOB_ID.out"
echo "  Job errors: logs/duke_job_$JOB_ID.err"
echo "  Duke log: logs/<TIMESTAMP>/<TIMESTAMP>_duke_run.log"
echo "========================================"

# ============================================================================
# TROUBLESHOOTING:
# ============================================================================
#
# Job fails immediately:
#   cat logs/duke_job_<JOB_ID>.err
#   Check module loading and paths
#
# Job runs but Duke fails:
#   cat logs/<TIMESTAMP>/<TIMESTAMP>_duke_run.log
#   ls -lh result_duke/module_data/
#
# Out of memory:
#   Increase: -l mem=16G
#   Enable: remove_intermediate = TRUE in duke_run_myriad.R
#
# Out of temp space:
#   Increase: -l tmpfs=100G
#
# Resume from failure:
#   Set resume = TRUE in duke_run_myriad.R
#   Delete broken cache: rm result_duke/temp/*.RData
#   Resubmit: qsub duke_myriad_job.sh
#
# ============================================================================
