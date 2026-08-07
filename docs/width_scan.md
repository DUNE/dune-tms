# Width-scan validation

`scripts/Validation/width_scan.py` runs the same edep-sim input files through
`ConvertToTMSTree` for a list of containment settings.  It creates a separate
TOML and ROOT output for each point, then writes `summary.csv` and
`summary.json` with candidate, track, hit, truth-denominator, and X-view Hough
rejection-stage counts.

After setting up the normal DUNE-TMS/ROOT environment and building this branch,
a 5,000-event single-file scan can be launched with:

```bash
python3 scripts/Validation/width_scan.py \
  --input /path/to/input.root \
  --output-dir width_scan_5000 \
  --max-events 5000 \
  --suite pr289
```

The helper automatically detects either `bin/ConvertToTMSTree.exe` (the Make
build) or `bin/ConvertToTMSTree` (the CMake-style name). Use `--reco-exe` only
to select another binary.

`--suite pr289` runs five exact, global-config comparisons: `legacy_ab`
(the diagnostic branch's deliberate pre-#289 X-only compatibility path),
`corrected_nominal` (Hough 1.5, extrapolation 1), `hough_4x` (6, 1),
`extrapolation_4x` (1.5, 4), and `retuned_4x` (6, 4). It enables compact
`Hough_Diagnostics` output but not the much larger hit snapshots. This is the
useful minimum scan for locating whether the loss is at initial X seeding,
DBSCAN/A*, or later 3-D reconstruction.

For an extra point, use `--point NAME,HOUGH_HALF_WIDTH,EXTRAPOLATION_MULTIPLIER`
and add `,true` only to request the legacy X compatibility path. These settings
are intentionally shared between X and Y in the corrected reconstruction; the
old four per-view knobs no longer exist.

The most useful comparisons are:

- `X attempts`, `X accepted`, `X no seed`, `X no DBSCAN`, `X A* empty`, and
  `X too short`: the actual Hough failure point;
- `X lines` and `Y lines`: where a setting changes candidate formation;
- `tracks`, `track slices`, and X/Y track hits: final reconstruction yield and
  view balance;
- truth primary muons touching the TMS: a stable denominator for a fixed input.

`tracks / truth #mu` is deliberately labeled as a relative scan diagnostic,
not a one-to-one matched efficiency.  A physics efficiency requires an explicit
reco-to-truth matching definition, which this lightweight validation does not
invent.
