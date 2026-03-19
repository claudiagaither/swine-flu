#!/bin/bash
#SBATCH -J treeannotator
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=128g
#SBATCH -t 12:00:00
#SBATCH -o treeannotator_%j.out
#SBATCH -e treeannotator_%j.err

module purge
module load beast/1.10.4

# Give Java most of the allocated memory
export _JAVA_OPTIONS="-Xmx120g -XX:+UseG1GC"

# Paths — update these to match your files
INTREES="combined_epifluHA_2010.trees"
OUTTREE="mcc_combined_epifluHA_2010.trees"

# Burn-in: 4491 trees (= 44906000 states from your run)
treeannotator \
    -burnin 4491 \
    -heights ca \
    "$INTREES" \
    "$OUTTREE"
