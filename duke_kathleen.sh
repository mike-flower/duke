#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_01_160c_2G_24h
#$ -l h_rt=24:00:00
#$ -pe mpi 160
#$ -l mem=2G
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ==============================================================================
# Duke Pipeline - Old Kathleen Job Script (CLI Version)
# ==============================================================================
# 
# This script runs Duke using the command-line interface (duke_cli.R)
# Edit the parameters in the ./duke command below
#
# Cluster: Old Kathleen (UCL) - Being phased out
# Scheduler: Grid Engine (SGE)
# Parallel environment: mpi (distributed)
# Cores: 80 (2 nodes × 40 cores)
# Memory: 4GB per core = 320GB total
# Runtime: 48 hours
#
# CRITICAL NOTES:
# - Old Kathleen requires -pe mpi (not -pe smp)
# - Request cores in multiples of 40 (node size)
# - No tmpfs support (tmpdir not available)
# - Being replaced by New Kathleen (Slurm)
#
# Expected performance (80 cores):
#   ~1000 samples: 15-18 hours
#   ~2000 samples: 30-36 hours
#
# For memory issues with very large datasets, see:
#   duke_kathleen_160c_80t_cli.sh (requests 160 cores, uses 80 threads)
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
  --threads 80 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE

# ==============================================================================
# COMMON CUSTOMIZATIONS
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

# Run specific modules only:
# ./duke --run_modules 5,6,7 --resume TRUE ...

# Custom repeat parameters:
# ./duke --rpt_pattern CTG --rpt_max_mismatch 1 ...

# ==============================================================================
# MEMORY TROUBLESHOOTING
# ==============================================================================

# If you encounter "Cannot allocate memory" errors with large datasets,
# use the memory-optimized version instead:
#   qsub duke_kathleen_160c_80t_cli.sh
#
# That script requests 160 cores but uses only 80 threads, providing
# a large memory buffer (320GB / 80 threads = 4GB per thread)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

echo ""
echo "Duke pipeline complete!"
echo "Results in: ~/Scratch/results/my_results"
