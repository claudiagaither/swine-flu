#!/bin/bash

#SBATCH --job-name=beast_run
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err
#SBATCH --time=168:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=franc@email.unc.edu

module load beast/1.10.4

# Save state every 10M steps so runs can be resumed with -resume
beast -threads $SLURM_CPUS_PER_TASK \
      -save_state epiflu_HA2010_v6.state \
      -save_every 10000000 \
      epiflu_HA2010_v6.xml
