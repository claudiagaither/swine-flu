#!/bin/bash
# ============================================================================
#  Resume incomplete BEAST runs — US_1990_v2 (12 trees)
# ============================================================================
#  Each tree folder has its own XML and .state file.
#  Uses -load_state to resume from the saved checkpoint.
#
#  Usage: bash resume_beast.sh
# ============================================================================

BASEDIR="/work/users/f/r/franc/swine_flu/BEAST_runs/US_1990_v2"
WALLTIME="72:00:00"

# Runs to resume (comment out or remove any that are already done)
RUNS=(01 02 03 04 05 06 07 08 09 10 11 12)

for RUN in "${RUNS[@]}"; do
    RUNDIR="${BASEDIR}/tree${RUN}"

    # Auto-detect the XML file
    XML=$(ls "${RUNDIR}"/*.xml 2>/dev/null | head -n 1)
    if [ -z "$XML" ]; then
        echo "WARNING: No .xml found in ${RUNDIR}, skipping tree${RUN}"
        continue
    fi
    XMLBASE=$(basename "$XML")

    # Auto-detect the .state file
    STATE=$(ls "${RUNDIR}"/*.state 2>/dev/null | head -n 1)
    if [ -z "$STATE" ]; then
        echo "WARNING: No .state file found in ${RUNDIR}, skipping tree${RUN}"
        continue
    fi
    STATEBASE=$(basename "$STATE")

    echo "Submitting resume for tree${RUN} (xml: ${XMLBASE}, state: ${STATEBASE})..."

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=beast_resume_t${RUN}
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
beast -load_state ${STATEBASE} -force_resume -overwrite -threads \$SLURM_CPUS_PER_TASK ${XMLBASE}
EOF

done

echo ""
echo "All jobs submitted. Check with: squeue -u franc"
