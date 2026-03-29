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
INTREES="/work/users/f/r/franc/swine_flu/combine_trees/H3N2_2010_v6/comb_epiflu_HA2010_v6.trees"
OUTTREE="mcc_comb_epiflu_HA2010_v6.trees"

# Burn-in: 5917 trees (= 5917000 states from your run)
treeannotator \
    -burnin 5917 \
    -heights ca \
    "$INTREES" \
    "$OUTTREE"
