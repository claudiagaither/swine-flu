#!/bin/bash

#SBATCH --job-name=flu_mafft_align
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=flu_align_%j.out
#SBATCH --error=flu_align_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@unc.edu   # <-- update this

# ── Load MAFFT ──────────────────────────────────────────────────────────────
module purge
module load mafft

# ── Paths ────────────────────────────────────────────────────────────────────
# Update INPUT to wherever you put your FASTA on Longleaf (e.g. /work/users/o/n/onyen/)
INPUT="/work/users/franc/swine_flu/gisaid_epiflu_sequence.fasta"
OUTPUT="/work/users/franc/swine_flu/aligned_mafft.fasta"

# ── Deduplicate sequence names before aligning ───────────────────────────────
# (mirrors the make.unique() step in your R script)
DEDUPED="${INPUT%.fasta}_deduped.fasta"

python3 - <<PYEOF
import sys
from collections import defaultdict

input_file  = "$INPUT"
output_file = "$DEDUPED"

counts = defaultdict(int)
records = []

with open(input_file) as fh:
    header, seq = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))

with open(output_file, "w") as out:
    for name, sequence in records:
        counts[name] += 1
        unique_name = name if counts[name] == 1 else f"{name}.{counts[name]-1}"
        out.write(f">{unique_name}\n{sequence}\n")

print(f"Written {len(records)} sequences to {output_file}")
PYEOF

# Stop immediately if the deduplication step failed
if [ $? -ne 0 ] || [ ! -s "$DEDUPED" ]; then
    echo "ERROR: Python deduplication failed or produced an empty file. Check paths above."
    exit 1
fi

# ── Run MAFFT ────────────────────────────────────────────────────────────────
# --auto   : lets MAFFT choose the best strategy for your data size
# --thread : use all allocated CPUs
# For very large datasets (>10k seqs) consider --retree 1 or --parttree
mafft \
    --auto \
    --thread ${SLURM_CPUS_PER_TASK} \
    --reorder \
    "$DEDUPED" \
    > "$OUTPUT"

echo "Alignment complete: $OUTPUT"
