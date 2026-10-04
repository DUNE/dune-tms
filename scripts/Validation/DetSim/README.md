# Detector-simulation validation

`detsim_validation.py` compares two or more `ConvertToTMSTree` runs over the same input files at the hit level (no
tracking or slicing) and writes `report.md`, `metrics.json` and six plots into an output directory. The first run is
the reference.

```
python3 detsim_validation.py --out <dir> [--files 1-15] off=<run dir> on=<run dir>
```

Each run directory holds `<NNNNNNN>_RecoCandidates.root` and `<NNNNNNN>_Readout.root` per file number. Needs uproot,
awkward, numpy and matplotlib (at FNAL: the CVMFS python v3_9_15 setup plus a user install of those packages).

What it measures (muon-dominated hits: `TrueLeptonicEnergy / TrueHitE > 0.95`):

| Section | Quantity | Plot |
|---|---|---|
| 1 | readout hits and PE above threshold per spill; PE after fibers per MeV for muon and hadronic hits | |
| 2 | hit survival vs true path length (`TrueHitDx`) and vs true deposited energy | `survival.png` |
| 3 | reco PE mean, relative spread and fraction below 3 PE per band; integer-PE fraction | `pe_partial_fill.png` |
| 4 | mean reco - true hit time vs number of Geant4 contributions, full crossings | `time_bias_vs_contributions.png` |
| 5 | one particle with more than one final hit in a bar and slice, by time gap (under the 120 ns window = merge failure) | `merge_gaps.png` |
| 6 | plane coverage: fraction of (particle, plane) deposits with a surviving hit, by deposited energy | `plane_coverage.png` |
| 7 | light provenance (`UseResponseElements = true` only): photon count vs PE, most-light vs most-energy particle | `light_share.png` |

`TrueHitDx` is not written by default; without it, sections 2-4 use deposited-energy bands instead of path length.

## Stage-by-stage scans (`stages/`)

Controlled checks of each simulation stage on its own, legacy vs new pipeline: synthetic crossings of one known bar
through the real `TMS_Event::FinalizeEvent()` (`ArtificialResegmentationTest` and `DetSimStageScan`), reduced to a
self-contained HTML page.

```
stages/run_stage_scans.sh <build dir> <edep-sim file for the geometry> <out>/scans      # ~1 min on 16 cores
python3 stages/stage_page_data.py --scans <out>/scans [--spills <new-pipeline run> --files 1-20] \
    [--metrics <detsim_validation.py metrics.json>] --out <out>/page_data.json
python3 stages/build_stage_page.py <out>/page_data.json notes.json <out>/detsim_stages.html
```

`notes.json` holds the page's prose (keys INTRO, META, N1-N6, N4X, OPEN, REGEN); it quotes numbers, so check it against
a new `page_data.json`.
