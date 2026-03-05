#!/bin/bash

#SBATCH --job-name=beast_run
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err
#SBATCH --time=96:00:00          # adjust as needed,  BEAST runs can be long
#SBATCH --mem=16G                # adjust based on your dataset size
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# --- Check/load BEAST module--- 
module load beast/1.10.4  # you may need to adjust the module name/version

# ---  Run BEAST--- 
beast -threads $SLURM_CPUS_PER_TASK epiflu_HA_2010.xml
