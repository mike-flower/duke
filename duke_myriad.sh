#!/bin/bash -l
#$ -N duke_pipeline
#$ -wd /home/skgtmdf/Scratch/bin/duke
#$ -o logs/duke_job_$JOB_ID.out
#$ -e logs/duke_job_$JOB_ID.err
#$ -l h_rt=48:00:00
#$ -pe smp 4
#$ -l mem=8G
#$ -l tmpfs=50G

# Load ALL required modules
module purge
module load r/recommended                   # R 4.2.0
module load samtools/1.11/gnu-4.9.2         # Samtools 1.11

# Set environment variables
export R_LIBS_USER=~/R/library
export PATH=$HOME/Scratch/bin/minimap2:$PATH

# Run Duke
cd ~/Scratch/bin/duke
Rscript duke_run_myriad.R
