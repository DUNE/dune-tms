#!/bin/bash
# Regenerates a "before this round's geometry fixes" and "after" ConvertToTMSTree
# output from the SAME input spill file, to check whether TMS_Bar's BarNumber
# formula change (TMS_Const::TMS_Start_Exact -> TMS_Geom survey-derived start)
# altered output values. Run this INSIDE the Apptainer/SL7 build container.
#
# Deliberately stashes only the five source files this round's geometry fixes
# touched -- NOT CMakeLists.txt/setup.sh, which have their own unrelated local
# changes that the build here depends on.
#
# Usage: from the repo root:
#   ./scripts/Validation/regen_before_after.sh /path/to/input.EDEPSIM_SPILLS.root
#
# IMPORTANT: ConvertToTMSTree's output-filename CLI argument is cosmetic only --
# TMS_TreeWriter/TMS_ReadoutTreeWriter always derive the real output name from the
# INPUT filename (see TMS_TreeWriter.cpp:17-33), unless overridden via the
# ND_PRODUCTION_TMSRECO_OUTFILE / ND_PRODUCTION_TMSRECOREADOUT_OUTFILE env vars.
# This script sets those explicitly for before/after so the two runs don't
# silently overwrite each other's output on the same auto-derived filename.
#
# Produces ./before_reco.root and ./after_reco.root (the file with Branch_Lines/
# Reco_Tree/Truth_Info -- what compare_bar_numbers.C reads), plus
# ./before_readout.root / ./after_readout.root (unused by the comparison, but
# still separated so they don't collide either).

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 /path/to/input.EDEPSIM_SPILLS.root"
  exit 1
fi
INPUT="$1"
if [ $# -gt 1 ]; then
  echo "Note: ignoring extra argument(s) '${*:2}' -- this script controls output naming itself via ND_PRODUCTION_TMSRECO_OUTFILE."
fi
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

GEOM_FILES=(
  src/TMS_Bar.cpp
  src/TMS_EventViewer.cpp
  src/TMS_Geom.cpp
  src/TMS_Kalman.cpp
  src/TMS_Reco.cpp
)

# Refuse to run with a dirty stash already present, to avoid confusing pop order
if git stash list | grep -q "regen_before_after.sh"; then
  echo "A stash from a previous run of this script is still present. Resolve it first:"
  git stash list
  exit 1
fi

echo "=== Stashing this round's geometry fixes (5 files) ==="
git stash push -m "regen_before_after.sh: geometry fixes" -- "${GEOM_FILES[@]}"

echo "=== Building BEFORE (pre-fix) ==="
cmake --build build -- -j"$(nproc)"

echo "=== Running BEFORE (pre-fix) ==="
ND_PRODUCTION_TMSRECO_OUTFILE="$REPO_ROOT/before_reco.root" \
  ND_PRODUCTION_TMSRECOREADOUT_OUTFILE="$REPO_ROOT/before_readout.root" \
  build/app/ConvertToTMSTree "$INPUT"

echo "=== Restoring this round's geometry fixes ==="
git stash pop

echo "=== Building AFTER (post-fix) ==="
cmake --build build -- -j"$(nproc)"

echo "=== Running AFTER (post-fix) ==="
ND_PRODUCTION_TMSRECO_OUTFILE="$REPO_ROOT/after_reco.root" \
  ND_PRODUCTION_TMSRECOREADOUT_OUTFILE="$REPO_ROOT/after_readout.root" \
  build/app/ConvertToTMSTree "$INPUT"

echo "=== Comparing bar-number branches ==="
root -b -q "scripts/Validation/compare_bar_numbers.C(\"before_reco.root\", \"after_reco.root\")"

echo
echo "Done. Now run:"
echo "  diff before_scan/RecoHitBar.txt after_scan/RecoHitBar.txt"
echo "  diff before_scan/RecoTrackKalmanBarView.txt after_scan/RecoTrackKalmanBarView.txt"
echo "  diff before_scan/TrueHitBar.txt after_scan/TrueHitBar.txt"
echo "No output from any of the three diffs means that branch is byte-identical before/after."
