#!/bin/bash
#SBATCH --job-name=logcombiner
#SBATCH --output=logcombiner_%j.out
#SBATCH --error=logcombiner_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# ============================================================================
#  Combine 8 BEAST .trees files with per-file burn-in using LogCombiner
# ============================================================================
#  Usage:
#    1. Edit the TREES array with your .trees file paths
#    2. Edit the BURNINS array with burn-in states for each file (check Tracer)
#    3. Edit OUTPUT_FILE with your desired output filename
#    4. Submit: sbatch combine_trees.sh
# ============================================================================

# --- Load BEAST (adjust module name if needed) ---
module purge
module load beast/1.10.4
# If your Longleaf module name differs, try:
#   module load beast/whatever version

# --- Configuration ---

# Input .trees files (edit these paths)
TREES=(
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree1/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree2/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree3/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree4/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree5/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree6/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree7/epiflu_HA2010_v6.trees.txt"
    "/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6/tree8/epiflu_HA2010_v6.trees.txt"
)

# Burn-in states for each file — SET THESE FROM TRACER
# (must be in the same order as TREES above)
BURNINS=(
    50000000   # run1 burn-in
    50000000   # run2 burn-in
    50000000   # run3 burn-in
    50000000   # run4 burn-in
    50000000   # run5 burn-in
    50000000   # run6 burn-in
    50000000   # run7 burn-in
    50000000   # run8 burn-in
)

# Resample every N states (20,000 = half the frequency of original 10,000)
RESAMPLE=20000

# Final combined output
OUTPUT_FILE="comb_epiflu_HA2010_v6.trees"

# Temp directory for intermediate per-file outputs
TMPDIR_WORK=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/logcombiner_XXXXXX")

# ============================================================================

echo "============================================"
echo "  LogCombiner — $(date)"
echo "  Working temp dir: ${TMPDIR_WORK}"
echo "============================================"

# --- Sanity checks ---
if [ ${#TREES[@]} -ne ${#BURNINS[@]} ]; then
    echo "ERROR: Number of tree files (${#TREES[@]}) does not match number of burn-in values (${#BURNINS[@]})."
    exit 1
fi

for f in "${TREES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Tree file not found: $f"
        exit 1
    fi
done

# --- Step 1: Strip per-file burn-in and resample each file individually ---
PROCESSED=()

for i in "${!TREES[@]}"; do
    INFILE="${TREES[$i]}"
    BURNIN="${BURNINS[$i]}"
    OUTFILE="${TMPDIR_WORK}/run$((i+1))_processed.trees"

    echo ""
    echo "Processing [$((i+1))/${#TREES[@]}]: ${INFILE}"
    echo "  Burn-in: ${BURNIN} states | Resample: ${RESAMPLE}"

    logcombiner -trees \
        -burnin "${BURNIN}" \
        -resample "${RESAMPLE}" \
        "${INFILE}" \
        "${OUTFILE}"

    if [ $? -ne 0 ]; then
        echo "ERROR: LogCombiner failed on ${INFILE}"
        rm -rf "${TMPDIR_WORK}"
        exit 1
    fi

    PROCESSED+=("${OUTFILE}")
    echo "  -> ${OUTFILE}"
done

# --- Step 2: Combine all processed files (burn-in already removed) ---
echo ""
echo "============================================"
echo "  Combining ${#PROCESSED[@]} processed files..."
echo "============================================"

logcombiner -trees \
    -burnin 0 \
    "${PROCESSED[@]}" \
    "${OUTPUT_FILE}"

if [ $? -ne 0 ]; then
    echo "ERROR: Final LogCombiner combine step failed."
    rm -rf "${TMPDIR_WORK}"
    exit 1
fi

# --- Cleanup ---
rm -rf "${TMPDIR_WORK}"

echo ""
echo "============================================"
echo "  Done! Output: ${OUTPUT_FILE}"
echo "  Finished at: $(date)"
echo "============================================"
