#!/bin/bash
# ============================================================================
#  Resume incomplete BEAST runs
# ============================================================================
#  This script submits a separate SLURM job for each run that needs resuming.
#  BEAST's -resume flag picks up from where each run left off.
#
#  Current status (from logs):
#    Run 1:  91.5M    Run 5: 100.0M (done)
#    Run 2:  93.7M    Run 6:  92.3M
#    Run 3:  91.4M    Run 7:  73.4M  <- needs extra time
#    Run 4:  98.6M    Run 8:  35.5M  <- needs extra time
#
#  Usage: bash resume_beast.sh
# ============================================================================

BASEDIR="/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6"
XML="epiflu_HA2010_v6.xml"

# Runs to resume (skip run 5 — already at 100M)
RUNS=(1 2 3 4 6 7 8)

for RUN in "${RUNS[@]}"; do
    RUNDIR="${BASEDIR}/tree${RUN}"

    # Sanity check
    if [ ! -f "${RUNDIR}/${XML}" ]; then
        echo "WARNING: XML not found in ${RUNDIR}, skipping run ${RUN}"
        continue
    fi

    # Runs 7 and 8 are further behind — give them the full 7 days
    if [ "$RUN" -eq 7 ] || [ "$RUN" -eq 8 ]; then
        WALLTIME="168:00:00"
    else
        WALLTIME="72:00:00"
    fi

    echo "Submitting resume for run ${RUN} (walltime: ${WALLTIME})..."

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=beast_resume_r${RUN}
#SBATCH --output=beast_resume_%j.out
#SBATCH --error=beast_resume_%j.err
#SBATCH --time=${WALLTIME}
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=franc@email.unc.edu

cd ${RUNDIR}
module load beast/1.10.4
beast -resume -threads \$SLURM_CPUS_PER_TASK ${XML}
EOF

done

echo ""
echo "All jobs submitted. Check with: squeue -u franc"
