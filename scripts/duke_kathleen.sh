#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_160c_2G_48h
#$ -l h_rt=48:00:00
#$ -pe mpi 160
#$ -l mem=2G
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_$JOB_ID.out
#$ -e logs/duke_$JOB_ID.err
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ==============================================================================
# Duke Pipeline - Old Kathleen Job Script
# Config: 160 cores requested, 80 threads used (2G/core = 320GB total)
# See README.md for resource guidance and Kathleen-specific notes
# ==============================================================================

# Environment
module purge
module load r/recommended
module load samtools/1.11/gnu-4.9.2
export R_LIBS_USER=~/R/library
export PATH=$HOME/Scratch/bin/minimap2:$PATH
mkdir -p logs

# Version check
echo "R:        $(R --version | head -1)"
echo "samtools: $(samtools --version | head -1)"
echo "minimap2: $(minimap2 --version)"
echo ""

# ==============================================================================
# EDIT PATHS BELOW
# ==============================================================================

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
