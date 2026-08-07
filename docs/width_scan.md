# Width-scan validation

`scripts/Validation/width_scan.py` runs the same edep-sim input files through
`ConvertToTMSTree` for a list of containment settings.  It creates a separate
TOML and ROOT output for each point, then writes `summary.csv` and
`summary.json` with candidate, track, hit, truth-denominator, and X-view Hough
rejection-stage counts.

The branch also writes compact truth-energy summaries for every X/Y 2-D Hough
candidate. The scan reports per-candidate mean completeness (the candidate
primary's energy divided by that particle's total cleaned energy in the same
view/slice), mean cleanliness (that primary energy divided by all truth energy
in the candidate), and primary-particle multiplicity. The latter is keyed by
run/spill/event/slice plus truth vertex and track ID, so a primary appearing in
two candidates is counted as a genuine duplicate reconstruction rather than a
cross-event coincidence.

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

The stage tables use a physics-facing denominator: primary muons whose actual
`Truth_Spill.BirthPosition` is inside the configured ND-LAr fiducial box and
that touch the TMS fiducial volume.  The scan deliberately does not use the
`LArFiducialStart` output branch for this selection.  A muon is reconstructed
at a stage only when it is the leading truth-energy
contributor to at least one candidate at that stage; the denominator is counted
once per spill, not once per time slice.  Completeness and cleanliness are
averaged over those target-muon candidates.

Truth summaries de-duplicate repeated copies of an identical reconstructed hit
(view, plane, bar, energy, time) before adding its truth-energy shares.  This
keeps completeness and cleanliness bounded by one when a reconstruction
working vector has carried the same hit more than once.
