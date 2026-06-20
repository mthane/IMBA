#!/usr/bin/env bash
# Same as run_python_tracker_example.ps1 — run from repo root after: pip install -e .
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
VIDEO="${1:-$REPO/legacy/data/example_video.mp4}"
if [[ ! -f "$VIDEO" ]]; then
  echo "Video not found: $VIDEO"
  echo "Copy example_video.mp4 into legacy/data/ or pass the path as the first argument."
  exit 1
fi
cd "$REPO"
python -m imba_tracker.cli --experiment-output --output-dir "$(dirname "$VIDEO")" "$VIDEO"
