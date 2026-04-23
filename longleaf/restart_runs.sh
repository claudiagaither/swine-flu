#!/bin/bash
# ============================================================================
#  Restart BEAST runs 7 and 8 from scratch (with state saving enabled)
# ============================================================================
#  Run 7 was at 73M, run 8 at 35M — restarting to get full 100M chains.
#  Old output files are backed up before restarting.
#
#  Usage: bash restart_runs_7_8.sh
# ============================================================================

BASEDIR="/nas/longleaf/home/franc/swine_flu/BEAST_runs/H3N2_2010_v6"
XML="epiflu_HA2010_v6.xml"

for RUN in 7 8; do
    RUNDIR="${BASEDIR}/tree${RUN}"

    if [ ! -f "${RUNDIR}/${XML}" ]; then
        echo "WARNING: XML not found in ${RUNDIR}, skipping run ${RUN}"
        continue
    fi

    # Back up old output files
    BACKUP="${RUNDIR}/old_run_backup"
    mkdir -p "${BACKUP}"
    echo "Backing up old output for run ${RUN} to ${BACKUP}/"
    mv "${RUNDIR}"/epiflu_HA2010_v6.trees.txt "${BACKUP}/" 2>/dev/null
    mv "${RUNDIR}"/epiflu_HA2010_v6.log.txt "${BACKUP}/" 2>/dev/null
    mv "${RUNDIR}"/beast_*.out "${BACKUP}/" 2>/dev/null
    mv "${RUNDIR}"/beast_*.err "${BACKUP}/" 2>/dev/null

    echo "Submitting fresh run ${RUN}..."

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=beast_r${RUN}
#SBATCH --output=beast_%j.out
#SBATCH --error=beast_%j.err
#SBATCH --time=168:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=franc@email.unc.edu

cd ${RUNDIR}
module load beast/1.10.4

beast -threads \$SLURM_CPUS_PER_TASK \
      -save_state epiflu_HA2010_v6.state \
      -save_every 10000000 \
      ${XML}
EOF

done

echo ""
echo "Jobs submitted. Check with: squeue -u franc"
