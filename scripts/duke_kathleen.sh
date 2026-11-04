#!/bin/bash -l
#$ -S /bin/bash
#$ -N duke_40c_4G_48h
#$ -l h_rt=48:00:00
#$ -pe mpi 40
#$ -l mem=4G
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_$JOB_ID.out
#$ -e logs/duke_$JOB_ID.err
#$ -M michael.flower@ucl.ac.uk
#$ -m bea

# ==============================================================================
# Duke - Kathleen job script (SGE)
# Duke runs samples in parallel on one node only, so request one node's cores: 40 cores,
# 4 GB per core, 48 h (check the node size on Kathleen and match it). Threads are set
# from $NSLOTS. Before the first run, edit -N, -wd, -M, the minimap2 PATH line and the
# paths below. Myriad (duke_myriad.sh) suits most runs.
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
