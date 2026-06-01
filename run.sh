#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Use the project virtualenv so python and memote resolve to .venv
if [ -f .venv/bin/activate ]; then
  source .venv/bin/activate
fi

echo "=== Running gem_annotate ==="
python -m scripts.gem_annotate

echo "=== Running memote ==="
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT="results/${TIMESTAMP}.html"
memote report snapshot --filename "$OUTPUT" model.xml

echo "=== Done: $OUTPUT ==="
