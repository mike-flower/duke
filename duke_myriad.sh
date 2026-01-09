#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_01_12c_8G_24h
#$ -l h_rt=24:00:00
#$ -pe smp 12
#$ -l mem=8G
#$ -l tmpfs=100G
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_$JOB_ID.out
#$ -e logs/duke_$JOB_ID.err
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ==============================================================================
# Duke Pipeline - Myriad Job Script (CLI Version)
# ==============================================================================
# 
# This script runs Duke using the command-line interface (duke_cli.R)
# Edit the parameters in the ./duke command below
#
# Cluster: Myriad (UCL)
# Scheduler: Grid Engine (SGE)
# Parallel environment: smp (shared memory)
# Cores: 36 (maximum for Myriad)
# Memory: 4GB per core = 144GB total
# Runtime: 48 hours
#
# Expected performance (36 cores):
#   ~100 samples: 2-3 hours
#   ~300 samples: 8-12 hours  
#   ~500 samples: 15-20 hours
#
# For larger datasets (1000+ samples), use Kathleen with 80+ cores
#
# ==============================================================================
# SUBMISSION INSTRUCTIONS:
# ==============================================================================
#
# 1. Make this script executable (one-time):
#    chmod +x duke_myriad.sh
#
# 2. Edit parameters in the ./duke command below
#
# 3. Submit the job:
#    qsub duke_myriad.sh
#
# 4. Monitor the job:
#    qstat -u $USER
#
# 5. Cancel job if needed:
#    qdel <JOB_ID>
#
# ==============================================================================

# Load required modules
module purge
module load r/recommended
module load samtools/1.11/gnu-4.9.2

# Set R library path
export R_LIBS_USER=~/R/library

# Add minimap2 to PATH
export PATH=$HOME/Scratch/bin/minimap2:$PATH

# Create logs directory if it doesn't exist
mkdir -p logs

# Verify modules loaded
echo "=== Loaded Modules ==="
module list
echo ""
echo "R version:"
R --version | head -1
echo ""
echo "samtools version:"
samtools --version | head -1
echo ""
echo "minimap2 version:"
minimap2 --version
echo ""

# Change to Duke directory
cd ~/Scratch/bin/duke

# ==============================================================================
# EDIT PARAMETERS BELOW
# ==============================================================================

# Run Duke with CLI
# Edit these paths for your analysis:

./duke \
  --dir_data /home/skgtmdf/Scratch/data/2025.12.17_pb_test/data \
  --dir_out /home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke \
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta \
  --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads 12 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE

# ==============================================================================
# COMMON CUSTOMISATIONS
# ==============================================================================

# Without adapter trimming:
# ./duke \
#   --dir_data ~/Scratch/data/my_experiment \
#   --dir_out ~/Scratch/results/my_results \
#   --path_ref ~/Scratch/refs/HTTset20/HTTset20.fasta \
#   --trim FALSE \
#   --threads 36

# Force re-run all modules:
# ./duke --resume FALSE ...

# Run specific modules only (e.g., re-plot after parameter change):
# ./duke --run_modules 5,6,7 --resume TRUE ...

# Downsample reads for testing:
# ./duke --downsample 1000 ...

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

echo ""
echo "Duke pipeline complete!"
echo "Results in: ~/Scratch/results/my_results"
