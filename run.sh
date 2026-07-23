#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Use the project virtualenv so python and memote resolve to .venv
if [ -f .venv/bin/activate ]; then
  source .venv/bin/activate
fi

echo "=== Running gem_annotate ==="
RESEARCH_ARGS=()
if [ -n "${IYALI26_RESEARCH_ROOT:-}" ]; then
  RESEARCH_ARGS=(--research-root "$IYALI26_RESEARCH_ROOT")
fi
python -m scripts.gem_annotate "${RESEARCH_ARGS[@]}"

# Sync the rebuilt model to the Genome-wide consumer repo (single source of truth
# = here). set -e above guarantees a failed rebuild aborts before this copy.
echo "=== Syncing model.xml to configured consumer ==="
CONSUMER_MODEL="${IYALI26_CONSUMER_MODEL:-}"
if [ -n "$CONSUMER_MODEL" ] && [ -d "$(dirname "$CONSUMER_MODEL")" ]; then
  cp model.xml "$CONSUMER_MODEL"
  echo "  synced -> $CONSUMER_MODEL  (md5 $(md5 -q model.xml))"
else
  echo "  IYALI26_CONSUMER_MODEL is unset or its directory is missing — skipping sync"
fi

echo "=== Running memote ==="
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT="results/${TIMESTAMP}.html"
memote report snapshot --solver glpk --filename "$OUTPUT" model.xml

echo "=== Done: $OUTPUT ==="
