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
CHUNK_COUNT="${IYALI26_DFBA_CHUNK_COUNT:-64}"
if [[ ! "$CHUNK_COUNT" =~ ^[1-9][0-9]*$ ]]; then
  echo "IYALI26_DFBA_CHUNK_COUNT must be a positive integer" >&2
  exit 2
fi
if (( SLURM_ARRAY_TASK_ID >= CHUNK_COUNT )); then
  echo "array task $SLURM_ARRAY_TASK_ID is outside chunk count $CHUNK_COUNT" >&2
  exit 2
fi
if [ "${IYALI26_DFBA_SOLVER:-glpk}" = "gurobi" ]; then
  module load gurobi/10.0.0
fi
cd "$IYALI26_REPO_ROOT"
TOL_ARGS=()
if [ -n "${IYALI26_DFBA_FEASIBILITY_TOL:-}" ]; then
  TOL_ARGS=(--feasibility-tol "$IYALI26_DFBA_FEASIBILITY_TOL")
fi
ALPHA_ARGS=()
if [ -n "${IYALI26_DFBA_ALPHA_REPLICATES:-}" ] && [ -z "${IYALI26_DFBA_ALPHA_SEED:-}" ]; then
  echo "IYALI26_DFBA_ALPHA_REPLICATES requires IYALI26_DFBA_ALPHA_SEED" >&2
  exit 2
fi
if [ -n "${IYALI26_DFBA_ALPHA_SEED:-}" ]; then
  ALPHA_ARGS=(
    --alpha-seed "$IYALI26_DFBA_ALPHA_SEED"
    --alpha-replicates "${IYALI26_DFBA_ALPHA_REPLICATES:-1}"
  )
fi
GENE_ARGS=()
if [[ ${IYALI26_DFBA_GENES+x} ]]; then
  if [[ -z "$IYALI26_DFBA_GENES" ]]; then
    echo "IYALI26_DFBA_GENES is set but empty" >&2
    exit 2
  fi
  IFS=',' read -r -a GENE_IDS <<< "${IYALI26_DFBA_GENES},"
  if (( ${#GENE_IDS[@]} == 0 )) || (( CHUNK_COUNT > ${#GENE_IDS[@]} )) || [[ -z "${GENE_IDS[0]}" ]]; then
    echo "targeted run needs 1 <= IYALI26_DFBA_CHUNK_COUNT <= number of IYALI26_DFBA_GENES" >&2
    exit 2
  fi
  for gene_id in "${GENE_IDS[@]}"; do
    if [[ -z "${gene_id//[[:space:]]/}" ]]; then
      echo "IYALI26_DFBA_GENES must not contain empty gene IDs" >&2
      exit 2
    fi
  done
  GENE_ARGS=(--genes "$IYALI26_DFBA_GENES")
fi
"$IYALI26_DFBA_PYTHON" -m scripts.gem_annotate.quinone_dfba_essentiality \
  --research-root "$IYALI26_RESEARCH_ROOT" --output-dir "$OUT" --run-id "$RUN_ID" \
  --solver "${IYALI26_DFBA_SOLVER:-glpk}" --uracil-mode "${IYALI26_DFBA_URACIL_MODE:-finite_batch}" \
  --dt "${IYALI26_DFBA_DT:-0.25}" "${TOL_ARGS[@]}" "${ALPHA_ARGS[@]}" "${GENE_ARGS[@]}" \
  --chunk-index "$SLURM_ARRAY_TASK_ID" --chunk-count "$CHUNK_COUNT"
