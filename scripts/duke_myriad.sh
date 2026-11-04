#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_12c_64G_48h
#$ -l h_rt=48:00:00
#$ -pe smp 12
#$ -l mem=64G
#$ -l tmpfs=500G
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_$JOB_ID.out
#$ -e logs/duke_$JOB_ID.err
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ==============================================================================
# Duke - Myriad job script (SGE)
# Resources: 12 cores, 64 GB per core (on Myriad, mem is per core: 768 GB in all), 48 h,
# 500 GB tmpfs. Duke runs on one node. A small run (tens of samples) needs far less,
# e.g. -pe smp 6 and -l mem=32G; a run that only rebuilds reports needs 1 core and 8G.
# Before the first run, edit: -N (job name), -wd (your Duke folder), -M (your email),
# the minimap2 PATH line, and the paths below. See README.md, "Running on Myriad".
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
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20_cagcaa.fasta \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads $NSLOTS \
  --resume TRUE

# Optional additions: add any of these on new lines after --resume TRUE, putting a
# backslash at the end of every line except the last.
#   Two tracts in one run:  --path_tracts /path/tracts.csv --composites "CAG_CCG_sum=CAG+CCG:sum"
#   UMIs:                   --umi TRUE   (flanks default to M13-F / ciosiF+1; see --help)
#   Allele structure:       --cluster TRUE
#   Adapter trimming:       --trim TRUE --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv
#   Check settings first:   --dry_run
