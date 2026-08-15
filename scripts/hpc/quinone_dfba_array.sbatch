#!/bin/bash -l
#SBATCH --job-name=iyali26-q9-dfba
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=0-63%16

set -euo pipefail
: "${IYALI26_REPO_ROOT:?set to the cloned repository}"
: "${IYALI26_RESEARCH_ROOT:?set to the external research workspace}"
: "${IYALI26_DFBA_PYTHON:?set to the tested Python executable on HPCC}"
RUN_ID="${IYALI26_DFBA_RUN_ID:?set one shared run ID before submission}"
OUT="$IYALI26_RESEARCH_ROOT/artifacts/results/quinone_dfba_essentiality/$RUN_ID"
if [ "${IYALI26_DFBA_SOLVER:-glpk}" = "gurobi" ]; then
  module load gurobi/10.0.0
fi
cd "$IYALI26_REPO_ROOT"
"$IYALI26_DFBA_PYTHON" -m scripts.gem_annotate.quinone_dfba_essentiality \
  --research-root "$IYALI26_RESEARCH_ROOT" --output-dir "$OUT" --run-id "$RUN_ID" \
  --solver "${IYALI26_DFBA_SOLVER:-glpk}" --chunk-index "$SLURM_ARRAY_TASK_ID" --chunk-count 64
