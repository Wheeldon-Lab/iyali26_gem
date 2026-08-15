#!/bin/bash -l
#SBATCH --job-name=iyali26-q9-dfba-merge
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail
: "${IYALI26_REPO_ROOT:?set to the cloned repository}"
: "${IYALI26_RESEARCH_ROOT:?set to the external research workspace}"
: "${IYALI26_DFBA_PYTHON:?set to the tested Python executable on HPCC}"
: "${IYALI26_DFBA_RUN_ID:?set one shared run ID before submission}"
cd "$IYALI26_REPO_ROOT"
"$IYALI26_DFBA_PYTHON" -m scripts.gem_annotate.quinone_dfba_essentiality --research-root "$IYALI26_RESEARCH_ROOT" \
  --output-dir "$IYALI26_RESEARCH_ROOT/artifacts/results/quinone_dfba_essentiality/$IYALI26_DFBA_RUN_ID" --merge
