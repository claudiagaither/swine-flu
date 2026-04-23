#!/bin/bash

#SBATCH --job-name=run_beast
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err
#SBATCH --time=240:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=franc@email.unc.edu

module load beast/1.10.4

# Save state every 100K steps so runs can be resumed with -resume
beast -threads $SLURM_CPUS_PER_TASK \
      -save_state 1990_2010.1like.state \
      -save_every 100000 \
      1990_2010.1like.xml
