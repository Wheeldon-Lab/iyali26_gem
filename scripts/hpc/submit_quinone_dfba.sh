#!/bin/bash
set -euo pipefail
: "${IYALI26_DFBA_RUN_ID:?export a shared run ID first}"
if [[ -n "${IYALI26_DFBA_ALPHAS:-}" && -n "${IYALI26_DFBA_ALPHA_SEED:-}" ]]; then
  echo "IYALI26_DFBA_ALPHAS and IYALI26_DFBA_ALPHA_SEED are mutually exclusive" >&2
  exit 2
fi
if [[ -z "${IYALI26_DFBA_ALPHAS:-}" && -z "${IYALI26_DFBA_ALPHA_SEED:-}" ]]; then
  IYALI26_DFBA_ALPHA_SEED="$(od -An -N8 -tu8 /dev/urandom | tr -d '[:space:]')"
  export IYALI26_DFBA_ALPHA_SEED
  echo "generated CoQ9 alpha seed: $IYALI26_DFBA_ALPHA_SEED"
fi
chunk_count="${IYALI26_DFBA_CHUNK_COUNT:-64}"
parallelism="${IYALI26_DFBA_ARRAY_PARALLELISM:-16}"
if [[ ! "$chunk_count" =~ ^[1-9][0-9]*$ ]] || [[ ! "$parallelism" =~ ^[1-9][0-9]*$ ]]; then
  echo "IYALI26_DFBA_CHUNK_COUNT and IYALI26_DFBA_ARRAY_PARALLELISM must be positive integers" >&2
  exit 2
fi
if [[ ${IYALI26_DFBA_GENES+x} ]]; then
  if [[ -z "$IYALI26_DFBA_GENES" ]]; then
    echo "IYALI26_DFBA_GENES is set but empty" >&2
    exit 2
  fi
  IFS=',' read -r -a genes <<< "${IYALI26_DFBA_GENES},"
  if (( ${#genes[@]} == 0 )) || [[ -z "${genes[0]}" ]]; then
    echo "IYALI26_DFBA_GENES must contain at least one comma-separated gene ID" >&2
    exit 2
  fi
  for gene_id in "${genes[@]}"; do
    if [[ -z "${gene_id//[[:space:]]/}" ]]; then
      echo "IYALI26_DFBA_GENES must not contain empty gene IDs" >&2
      exit 2
    fi
  done
  if (( chunk_count > ${#genes[@]} )); then
    chunk_count=${#genes[@]}
  fi
fi
export IYALI26_DFBA_CHUNK_COUNT="$chunk_count"
array_spec="0-$((chunk_count - 1))%$parallelism"
job_id=$(sbatch --parsable --array="$array_spec" scripts/hpc/quinone_dfba_array.sbatch)
sbatch --dependency="afterok:$job_id" scripts/hpc/quinone_dfba_merge.sbatch
echo "array job: $job_id"
