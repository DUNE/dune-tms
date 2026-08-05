# Width-scan validation

`scripts/Validation/width_scan.py` runs the same edep-sim input files through
`ConvertToTMSTree` for a list of transverse-containment settings.  It creates a
separate TOML and ROOT output for each point, then writes `summary.csv` and
`summary.json` with candidate, track, hit, and truth-denominator counts.

After setting up the normal DUNE-TMS/ROOT environment and building this branch,
a 5,000-event single-file scan can be launched with:

```bash
python3 scripts/Validation/width_scan.py \
  --input /path/to/input.root \
  --output-dir width_scan_5000 \
  --max-events 5000 \
  --point legacy_x,116.25,1.5,78,1 \
  --point widen_both_2x,3,3,2,2 \
  --point widen_both_4x,6,6,4,4
```

The automatic `default` point is `Hough X/Y = 1.5/1.5` and extrapolation X/Y
scales `1/1`.  `legacy_x` approximately recovers the pre-PR-289 X acceptance:
the Hough half-window needs about `1.5 * (bar length / bar width)`, while the
extrapolation multiplier needs about `bar length / bar width`.

The most useful comparisons are:

- `X lines`, `Y lines`, and `3D lines`: where a setting changes candidate
  formation;
- `tracks`, `track slices`, and X/Y track hits: final reconstruction yield and
  view balance;
- truth primary muons touching the TMS: a stable denominator for a fixed input.

`tracks / truth #mu` is deliberately labeled as a relative scan diagnostic,
not a one-to-one matched efficiency.  A physics efficiency requires an explicit
reco-to-truth matching definition, which this lightweight validation does not
invent.
