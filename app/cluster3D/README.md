# app/cluster3D — validation and diagnostic tools for the Cluster3D pipeline

These programs exercise the space-point reconstruction in `src/Cluster3D/`
(space points → DBSCAN + PCA → cluster linking → Kalman follower and split
step → shadow absorption and stitching; see `src/Cluster3D/README.md` for the
pipeline itself). None of them are part of production reconstruction, but
`Cluster3DRecoTruth` runs exactly the stage conversion runs
(`TMS_Cluster3DReco::Run()`); the older tools assemble their own variants of
it. Most read space points back out of a
`Reco_Tree` that `ConvertToTMSTree` already wrote and compare against the
truth stored in the same file.

The sources live here, but the binaries are built into `build/app/` like every
other app (see `app/CMakeLists.txt`), so existing scripts keep working.

**Geometry file.** Tools that fit tracks or use `TMS_Geom` need a file that
carries the `EDepSimGeometry` key (e.g. the input `*.EDEPSIM_SPILLS.root`).
The `*_RecoCandidates.root` output does not embed one.

**Environment hooks.** Configuration sweeps are done through environment
variables so settings can be compared without a rebuild. An unset variable
means the library default.

## Where to start

| If you want to… | Use |
|---|---|
| score the production Cluster3D reconstruction: every track and every true muon | `Cluster3DRecoTruth` |
| measure how often each true muon is found and how well it is fitted | `KalmanFollowerTruthEfficiency` |
| count *every* fitted track, including fakes and duplicates | `TrackFindingObjectTruth` |
| look at one slice in detail (every candidate, every χ²) | `KalmanFollowerSliceTest` |
| study the clustering stage alone | `ClusterTruthEfficiency` |

## Truth-validation tools

### `Cluster3DRecoTruth`
```
Cluster3DRecoTruth <geom.root> <reco.root> <tracks.csv> <muons.csv>
```
Object-first, hit-level validation of `TMS_Cluster3DReco::Run()` as
conversion runs it (the slice's hits for the hit-level fit and orphan pickup,
the X/Y time term), with no truth input. Truth only scores afterwards, from
the energy shares of the hits the tracks used:
- a track's owner is the particle with the largest summed share, if that is
  more than half the hits, otherwise the track is junk ("mixed");
- a true muon (\|PDG\| 13, at least 5 true hits in the slice) is found if it
  owns a track; its best track gives completeness and purity; further owned
  tracks are duplicates.

One row per track (owner, purity, hits per view and the owner's share of each,
start/end, stop reason, chi2, timing) and one row per slice × muon (outcome,
where its space points went among the DBSCAN clusters, and for a track ending
early, where the rest of the muon went). A muon cut across two slices has a
row in each, so count distinct muons by their main slice (most true hits).
- Reconstruction settings: `CLUSTER3D_LINK`, `CLUSTER3D_GRAPH` (and
  `_GRAPH_MIN_CLUSTER`, `_GRAPH_MIN_LAYERS`, `_GRAPH_MIN_HITS`),
  `CLUSTER3D_MIN_HITS`, `CLUSTER3D_SHADOW` (`_SHADOW_FRAC`, `_SHADOW_TOL`),
  `CLUSTER3D_STITCH` (`_STITCH_GAP`, `_STITCH_ANGLE`, `_STITCH_MISS_BASE`,
  `_STITCH_MISS_PER_M`), `CLUSTER3D_LINK_MAX_GAP`, `_LINK_MAX_ANGLE`,
  `_LINK_MISS_BASE`, `_LINK_MISS_PER_M`, and the follower's
  `CLUSTER3D_EXTEND` (`_EXTEND_GAP`, `_EXTEND_CHI2`),
  `CLUSTER3D_RESEED_FACTOR` (`_RESEED_SIGMA`, `_RESEED_PASSES`),
  `CLUSTER3D_RANGE_SEEDED`, `CLUSTER3D_RANGE_SEED_REL_SIGMA`,
  `CLUSTER3D_MAX_BEYOND_SEED`, `CLUSTER3D_STOP_ON_RANGE_OUT`,
  `CLUSTER3D_XYTIME_SIGMA`, `CLUSTER3D_XYTIME_GATE`, `CLUSTER3D_XYTIME_ALT`,
  `CLUSTER3D_RANGE_TABLES` (tracking), `CLUSTER3D_RANGE_TABLES_MOMENTUM` (range
  walk), `CLUSTER3D_RANGE_TABLES_SCALE` (stopping-power scale), `CLUSTER3D_EXPECTED_STOP`, `CLUSTER3D_MAX_GAP` (see
  `TMS_Cluster3DReco::Config` and `TMS_KalmanFollower::Config`). The
  conversion's own `[Recon.Cluster3D]` settings are not read: set
  `CLUSTER3D_LINK=1` to match a conversion with `LinkClusters = true`.
- `C3D_POINT_DUMP=1` also writes `<tracks.csv>.points.csv`: every track's
  chosen points with their transit-corrected X/Y time difference, both hits'
  truth, and whether a point built from two of the owner's hits was among the
  node's candidates (the view-mismatch study, 2026-09-28).

### `ClusterTruthEfficiency`
```
ClusterTruthEfficiency <geom.root> <reco.root> <muons.csv> <clusters.csv> [display_prefix]
```
Runs DBSCAN + PCA on every slice and matches track-like clusters to true muons
by plurality vote. Reports per-muon cluster efficiency, hit completeness and
purity. It defines the truth conventions every other tool here copies: loading
`Truth_Spill`, and collapsing a delta ray or other descendant to its muon
ancestor through the `Parent` chain.
- `CTE_CLUSTER_DETAIL_CSV=<path>` writes one row per cluster, track-like or
  not: PCA eigenvalues, per-layer transverse span, the top two truth owners and
  their mean times. Used for the merged-muon study
  (`reports/2026-09-24_caseH_timing_pca/`).

### `GraphTrackFinderTruthEfficiency`
```
GraphTrackFinderTruthEfficiency <geom.root> <reco.root> <muons.csv> [append 0|1]
```
Validates the three-stage finder (DBSCAN + PCA, then merged-cluster PCA, then
`TMS_GraphTrackFinder`) over the whole muon population. No Kalman fit.
- `LT_DEBUG_VGID=<vgid>` prints debug output for one muon.

### `KalmanFollowerTruthEfficiency`
```
KalmanFollowerTruthEfficiency <geom.root> <reco.root> <muons.csv> [append 0|1]
```
The main per-muon benchmark. It runs the same three stages, fits whichever
object is found with `TMS_KalmanFollower`, and scores the fit against truth:
completeness, purity, ambiguity-resolution accuracy and why the walk stopped.
The `kalman_strict_*` columns count a point as correct only when its X-hit and
Y-hit truth labels both belong to the target. The plain columns use one label
and so also credit ghosts, points pairing the target's hit with another
particle's.

Muon-first: it starts from each true muon, so it cannot see fake tracks. Use
`TrackFindingObjectTruth` for those.
- Follower settings: `KF_QP_REL_SIGMA`, `KF_RANGE_SEED`, `KF_STOP_ON_RANGEOUT`,
  `KF_MAX_HEAD_SKIP`, `KF_MAX_TRIPLETS`, `KF_RANK_BY_CONVERGENCE`,
  `KF_USE_TIME`, `KF_TIME_SIGMA`, `KF_TIME_GATE` (see
  `TMS_KalmanFollower::Config` for meanings and defaults).
- `KF_DUMP_HYPOTHESES=<path>` writes one row per seed hypothesis tried.
- `KF_DUMP_MISSED=<path>` writes one row per target plane the fit walked past
  without picking the right point.

### `TrackFindingObjectTruth`
```
TrackFindingObjectTruth <geom.root> <reco.root> <tracks.csv> <muons.csv>
```
Object-first validation: runs the reconstruction on every slice with no truth
input (DBSCAN, then `TMS_IterativeTrackFit` on every track-like cluster), then
classifies each fitted track against truth: clean muon, duplicate of an
already-found muon, owned by another particle, or made only of ghost points.
The muons CSV gives per-muon efficiency from the same output. Only track-like
clusters are fitted, so absolute efficiency is lower than the muon-first
tools report; the tool is meant for comparing variants and for fake rates.
- `ITF_SPLIT` = `0` off, `1` flagged clusters only (default), `2` every
  cluster.
- `ITF_MAX_TRACKS`, `ITF_MIN_HITS`, `ITF_MIN_SPLIT_HITS`, `ITF_RESTRICT_SPLIT`
  (see `TMS_IterativeTrackFit::Config`).
- `KF_MAX_HEAD_SKIP`, `KF_MAX_TRIPLETS`, `KF_USE_TIME`, `KF_TIME_SIGMA`,
  `KF_TIME_GATE`, as above.

## Single-case tools

### `KalmanFollowerSliceTest`
```
KalmanFollowerSliceTest <geom.root> <reco.root> [vgid] [trackid] [display.json]
```
Runs the full pipeline on the slice that holds the most of one true
particle's points, and prints the fit node by node with a truth check. Each
node is marked `correct`, `HALF` (a ghost with one coordinate from another
particle) or `WRONG`. The optional JSON feeds the event-display pages.
- `KF_DUMP_CANDIDATES=1` lists every candidate at every layer, with its χ²,
  time χ², X-hit and Y-hit truth labels, and time.
- `KF_DEBUG=1` prints the clusters that touch the target's own points.
- Same `KF_*` settings as `KalmanFollowerTruthEfficiency`.

### `GraphTrackFinderSliceTest`
```
GraphTrackFinderSliceTest <reco.root> [vgid] [trackid]
```
Runs `TMS_GraphTrackFinder` on one known-hard slice (by default the
shower-contaminated muon that motivated the graph search) and reports its
candidate paths against truth.

### `GraphTrackFinderSweep`
```
GraphTrackFinderSweep <targets.csv (file,vgid,trackid)> <results.csv>
```
Runs the graph search over a list of muons DBSCAN + PCA failed on, reporting
how many it recovers. Used to tune `TMS_GraphTrackFinder::Config`.

## Event-display and dump tools

### `RealEventClusterDisplay` / `SpillClusterDisplay`
```
RealEventClusterDisplay <geom.root> <reco.root> <entry> <max_dz_mm> <base_transverse_mm> <transverse_per_dz> <min_points> <out.csv>
SpillClusterDisplay     <geom.root> <reco.root> <spill> <max_dz_mm> <base_transverse_mm> <transverse_per_dz> <min_points> <out.csv>
```
The DBSCAN tolerance arguments are `TMS_SpacePointDBScan::Tolerance`
(defaults: 270, one bar pitch, 0.55).
Run the real DBSCAN + PCA on one slice (or on every slice of one spill,
clustered per slice as in reconstruction), and write per-point cluster labels
and per-cluster PCA to CSV for 3D display.

### `DumpTrueTrajectory`
```
DumpTrueTrajectory <geom.root> <out.json> <vertexId>:<trackId> [...]
```
Writes the Geant4 trajectory points of one or more true particles inside the
TMS, read from the raw `TG4Event`. It uses no reconstruction at all; this is
the "true trajectory" overlay in the Kalman event displays.

### `DumpSpillTrajectories`
```
DumpSpillTrajectories <spill_file.root> <spillNumber> <out.json>
```
Writes every charged Geant4 trajectory that reaches the TMS in one spill, with
its PDG code, parent, initial momentum and all of its points (times relative
to the spill start), read from the raw `TG4Event`. Interactions are grouped
into spills exactly as `ConvertToTMSTree` does, so `spillNumber` matches the
converter's `SpillNumber`. Built for a whole-spill truth display.

### `DumpHitPE`
```
DumpHitPE <input.EDEPSIM_SPILLS.root> <hits.csv>
```
Writes each hit's true deposited energy, simulated photoelectrons and
pedestal-suppression flag, taken just after the detector-response chain and
before slicing or reconstruction. Built for the pedestal-threshold study; kept
here with the other tools from the same branch.

## Unit-level tests

### `SpacePointCluster_test`
```
SpacePointCluster_test <max_dz_mm> <base_transverse_mm> <transverse_per_dz> <min_points>
```
Checks `TMS_KDTree` radius queries against brute force, and runs
`TMS_SpacePointDBScan` on synthetic data (two Gaussian blobs and two lines,
on a 40 mm synthetic plane grid: `120 40 1.0 5` is the old default).

### `GraphTrackFinder_test`
Checks `TMS_GraphTrackFinder` on synthetic tracks.
