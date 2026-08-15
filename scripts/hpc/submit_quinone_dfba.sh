#!/bin/bash
set -euo pipefail
: "${IYALI26_DFBA_RUN_ID:?export a shared run ID first}"
job_id=$(sbatch --parsable scripts/hpc/quinone_dfba_array.sbatch)
sbatch --dependency="afterok:$job_id" scripts/hpc/quinone_dfba_merge.sbatch
echo "array job: $job_id"
