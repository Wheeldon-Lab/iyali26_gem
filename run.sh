#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Use the project virtualenv so python and memote resolve to .venv
if [ -f .venv/bin/activate ]; then
  source .venv/bin/activate
fi

echo "=== Running gem_annotate ==="
python -m scripts.gem_annotate

# Sync the rebuilt model to the Genome-wide consumer repo (single source of truth
# = here). set -e above guarantees a failed rebuild aborts before this copy.
echo "=== Syncing model.xml to Genome-wide (consumer repo) ==="
GW_MODEL="/Users/david/Desktop/Lab/Ian wheeldon/code/Genome-wide/model.xml"
if [ -d "$(dirname "$GW_MODEL")" ]; then
  cp model.xml "$GW_MODEL"
  echo "  synced -> $GW_MODEL  (md5 $(md5 -q model.xml))"
else
  echo "  Genome-wide dir not found — skipping sync"
fi

echo "=== Running memote ==="
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT="results/${TIMESTAMP}.html"
memote report snapshot --solver glpk --filename "$OUTPUT" model.xml

echo "=== Done: $OUTPUT ==="
