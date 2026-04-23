O#!/bin/bash

#SBATCH --job-name=beast_run
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err
#SBATCH --time=96:00:00          # adjust as needed — BEAST runs can be long
#SBATCH --mem=16G                # adjust based on your dataset size
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# --- Check/load BEAST module ---
module load beast2  # you may need to adjust the module name/version

# --- Run BEAST ---
beast2 -threads $SLURM_CPUS_PER_TASK H3N2_50_B277.xml
