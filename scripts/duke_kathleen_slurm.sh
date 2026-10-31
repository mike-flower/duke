#!/bin/bash -l
#SBATCH --job-name=duke_80t_2G_48h
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
# Duke 3.2.0 - New Kathleen Job Script (Slurm)
# Config: 80 tasks across 2 nodes, 2GB/task = 160GB total, 48h
# Note: submit with sbatch (not qsub); no chmod needed
# See README.md for resource guidance and SGE vs Slurm command reference
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
echo "Nodes: $SLURM_JOB_NUM_NODES  Tasks: $SLURM_NTASKS  Job: $SLURM_JOB_ID"
echo ""

cd ~/Scratch/bin/duke

# ==============================================================================
# EDIT PATHS BELOW
# Note: --threads is set automatically from $SLURM_NTASKS
# ==============================================================================

./duke \
  --dir_data /home/skgtmdf/Scratch/data/2025.12.17_pb_test/data \
  --dir_out /home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke \
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20_cagcaa.fasta \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads $SLURM_NTASKS \
  --resume TRUE

# Optional additions: add any of these on new lines after --resume TRUE, putting a
# backslash at the end of every line except the last.
#   Two tracts in one run:  --path_tracts /path/tracts.csv --composites "CAG_CCG_sum=CAG+CCG:sum"
#   UMIs:                   --umi TRUE   (flanks default to M13-F / ciosiF+1; see --help)
#   Allele structure:       --cluster TRUE
#   Adapter trimming:       --trim TRUE --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv
#   Check settings first:   --dry_run

DUKE_EXIT_CODE=$?
echo ""
if [ $DUKE_EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully ($(date), duration: $(date -ud "@$SECONDS" +%T))"
else
    echo "Pipeline FAILED with exit code $DUKE_EXIT_CODE ($(date))"
fi
exit $DUKE_EXIT_CODE
