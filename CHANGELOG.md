# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to a modified version of [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Geometry releases will be tagged as `Descriptive_tag_v_X.Y.Z`.

## [1.2.0] - Unreleased

### Added

- Adds `.github/copilot-instructions.md` to guide Copilot reviews with repo-specific conventions (search before claiming something is unused, cite file/line, separate bugs from style). ([#339](https://github.com/DUNE/dune-tms/pull/339))
- Adds per-stage wall-clock timers around the space-point builder and Cluster3D reconstruction in the converter, and adds build artifact patterns to `.gitignore`. ([#337](https://github.com/DUNE/dune-tms/pull/337))
- Introduces the full Cluster3D pipeline (DBSCAN clustering, cluster linking, Kalman-follower with hit-time gating, range/curvature momentum) running alongside legacy tracking and writing to `Reco_Tree_C3D`/`Truth_Info_C3D`; disabled by default, with significant efficiency gains when enabled (93.1% vs 82.6% ND-LAr muons found). ([#336](https://github.com/DUNE/dune-tms/pull/336))
- Introduces `TMS_SpacePoint` to match X- and Y-view hits from the same plane into 3D space points needed by the Cluster3D reconstruction; no effect on any existing output branch. ([#334](https://github.com/DUNE/dune-tms/pull/334))
- First scaffold of the detector-response redesign: adds a disabled `Sim.DetSim.UseResponseElements` config flag and an `ArtificialResegmentationTest` benchmark that exposes hit-time drift as a function of hit-segment count, with no change to simulation output. ([#331](https://github.com/DUNE/dune-tms/pull/331))
- Replaces the single energy-window slicer with a view-aware `PerViewTimeSlicer` that slices X and Y hits separately then links them via light-travel-corrected hit-time coincidences, reducing muons split across views from 5.3% to 1.1%. ([#329](https://github.com/DUNE/dune-tms/pull/329))
- Adds order-independent `GetPrimaryIdByEnergy()`/`GetVertexGlobalIdByEnergy()` to `TMS_TrueHit`, selecting the highest-energy contributor instead of whichever was first-merged, and exposes them as new `TrueHitPrimaryId`/`TrueHitVertexId` branches in `Truth_Info` (purely additive). ([#321](https://github.com/DUNE/dune-tms/pull/321))


### Changed

- Ports Cluster3D validation support onto `main` (configurable tree names, `--cluster3d` flag, FNAL-directory guard) and adds momentum-resolution and ND-LAr muon efficiency/duplicate-rate plots to `Tracking_Validation`. ([#341](https://github.com/DUNE/dune-tms/pull/341))
- Sets `[Recon.Cluster3D] Enabled = true` as the default config so the Cluster3D and legacy reconstructors run in parallel and both write output trees on every conversion. ([#338](https://github.com/DUNE/dune-tms/pull/338))
- Raises the 3D track-matching coincidence window from 16 to 30 ns to recover muons whose X and Y view hit times differ due to fiber-travel delay in the X/Y geometry, gaining ~8 percentage points of ND-LAr muon efficiency. ([#333](https://github.com/DUNE/dune-tms/pull/333))
- Reduces per-event loop overhead in `ConvertToTMSTree` by ~65% through algorithmic improvements, with output that is bit-for-bit identical to the previous version. ([#332](https://github.com/DUNE/dune-tms/pull/332))


### Fixed

- Corrects the comment-only `PrintSummary()` usage example in `TMS_StageTimer.h`, which incorrectly showed no arguments instead of the required `(totalSeconds, nEntries)`. ([#340](https://github.com/DUNE/dune-tms/pull/340))
- Fixes an off-by-one in the `GetMaterials()` material-walk loop that accumulated extra path beyond the endpoint, causing step-size-dependent material integrals; reduces the legacy range-momentum bias from +1.4% to −3.9% and corrects Cluster3D's range momentum. ([#335](https://github.com/DUNE/dune-tms/pull/335))
- Fixes `MergeCoincidentHits()` and `SimulateDeadtime()` to distinguish split X-bar halves by bar number rather than just `(z, NotZ)`, preventing incorrect hit merging and deadtime cross-contamination between physically separate channels (~1% more tracks reconstructed). ([#328](https://github.com/DUNE/dune-tms/pull/328))
- Ensures `VisibleE` is properly accumulated in `TMS_TrueHit::MergeWith()`, preventing visible energy from being silently dropped when true hits are combined. ([#327](https://github.com/DUNE/dune-tms/pull/327))
- Fixes `TMS_TrueHit::MergeWith()` to recompute `PEAfterFibers` as the sum of short- and long-path PE components rather than leaving the value stale after a merge. ([#326](https://github.com/DUNE/dune-tms/pull/326))
- Defaults `TMS_Hit::RecoX/Y` and their uncertainties to `-999` and initialises `nVertices = -999` in the slice constructor, eliminating uninitialised-memory reads that were silently corrupting five `Truth_Info` branches on every spill entry. ([#325](https://github.com/DUNE/dune-tms/pull/325))
- Corrects the validation reader headers to declare `StartPos`/`EndPos` as `[nTracks][4]` (matching the tree's actual layout since May 2025) rather than `[MAX_TRACKS][3]`, fixing silently shifted position reads for all tracks after the first in a slice. ([#324](https://github.com/DUNE/dune-tms/pull/324))
- Resets the Y-candidate iterator to `begin()` after `erase()` whenever the list is non-empty, eliminating a dangling-iterator dereference that made 3D track finding nondeterministically dependent on heap contents. ([#323](https://github.com/DUNE/dune-tms/pull/323))
- Guards `RunKalman()` and `runRTSSmoother()` against tracks with fewer than two Kalman nodes, eliminating a deterministic heap corruption and `TMatrixT` abort on certain spill files. ([#322](https://github.com/DUNE/dune-tms/pull/322))


## [1.1.0] - Unreleased

## [MAIN] - Unreleased
- Changed MergeTracks default to True.
- Hits trailing the end of a Hough Track can now be added to the track they 'should' correspond to. ON by default, switchable in config.
- Prevent 3D track matching from looping indefinitely when the last hit in a required view lies outside the detector bounds.

### Added

- Added run-aware global vertex identifiers for truth tracking ([#247](https://github.com/DUNE/dune-tms/pull/247)).
- Restored the GitHub Actions build workflow ([#259](https://github.com/DUNE/dune-tms/pull/259)).
- Added `nTrueParticles` back to `Truth_Info` ([#282](https://github.com/DUNE/dune-tms/pull/282)).
- Added the missing Y-view line-candidate output and corrected X/Y branch names ([#284](https://github.com/DUNE/dune-tms/pull/284)).
- Added leading truth-contributor information for reconstructed hits ([#303](https://github.com/DUNE/dune-tms/pull/303)).
- Added a geometry survey that derives detector bounds, planes, modules, and configurable volume names from the loaded GDML ([#298](https://github.com/DUNE/dune-tms/pull/298)).
- Added a truth side table and truthless reconstructed-hit construction to separate simulation truth from reconstructed detector data ([#301](https://github.com/DUNE/dune-tms/pull/301)).
- Retained sampled sensor-photon arrival times, including their source hit and fiber-path metadata, as groundwork for future readout-window handling ([#318](https://github.com/DUNE/dune-tms/pull/318)).

### Changed

- Made DBSCAN reconstruction use its configured values ([#244](https://github.com/DUNE/dune-tms/pull/244)).
- Removed reflection from the default readout configuration ([#246](https://github.com/DUNE/dune-tms/pull/246)).
- Corrected plane numbering and updated dependent geometry behavior ([#249](https://github.com/DUNE/dune-tms/pull/249)).
- Made A* reconstruction use its configured values ([#268](https://github.com/DUNE/dune-tms/pull/268)).
- Retuned the default reconstruction configuration for efficiency ([#271](https://github.com/DUNE/dune-tms/pull/271)).
- Made the better Kalman chi-squared-per-degree-of-freedom fit canonical for derived track variables ([#269](https://github.com/DUNE/dune-tms/pull/269)).
- Built contiguous bar numbers from the detector geometry instead of coordinate assumptions ([#280](https://github.com/DUNE/dune-tms/pull/280)).
- Increased plane-matching limits to recover reconstruction efficiency ([#295](https://github.com/DUNE/dune-tms/pull/295)).
- Used the transverse X-bar width for Hough containment and made the widened window continuous ([#289](https://github.com/DUNE/dune-tms/pull/289)).
- Split detector simulation from signal processing to provide a reconstruction entry point for real detector data ([#300](https://github.com/DUNE/dune-tms/pull/300)).
- Updated the DUNE UPS product versions, build paths, and ROOT linking for the current build environment ([#302](https://github.com/DUNE/dune-tms/pull/302)).
- Aligned the CMake package version and ROOT output metadata version at `1.1.0`.
- Changed `MergeTracks` to default to true.
- Enabled adding trailing Hough-track hits to their corresponding track by default, with a configuration switch.
- Seeded detector-response randomness from the edep-sim run and event IDs, avoiding repeated streams across independent jobs and making an event's response stable despite earlier slice-count changes ([#319](https://github.com/DUNE/dune-tms/pull/319)).

### Fixed

- Used the configured magnetic-field boundary in the Kalman fit ([#253](https://github.com/DUNE/dune-tms/pull/253)).
- Handled time-slicer bounds correctly when time slicing is disabled ([#257](https://github.com/DUNE/dune-tms/pull/257)).
- Added a fallback when a reconstructed track length is unavailable ([#264](https://github.com/DUNE/dune-tms/pull/264)).
- Initialized and copied reconstructed track values consistently ([#262](https://github.com/DUNE/dune-tms/pull/262)).
- Cleared reconstructed output arrays consistently between events ([#263](https://github.com/DUNE/dune-tms/pull/263)).
- Preserved the invalid-energy sentinel in invalid truth momenta ([#260](https://github.com/DUNE/dune-tms/pull/260)).
- Stored truth primary flags with the correct type ([#275](https://github.com/DUNE/dune-tms/pull/275)).
- Wrote truth end positions through the common position writer ([#274](https://github.com/DUNE/dune-tms/pull/274)).
- Fixed signed/unsigned comparisons in track-direction calculations ([#283](https://github.com/DUNE/dune-tms/pull/283)).
- Avoided repeatedly fetching the same single-interaction gRoo entry ([#276](https://github.com/DUNE/dune-tms/pull/276)).
- Corrected LAr fiducial truth flags using the shared geometry predicates ([#297](https://github.com/DUNE/dune-tms/pull/297)).
- Removed an invalid deadtime-map access and the duplicate hit-merging workaround it had motivated ([#299](https://github.com/DUNE/dune-tms/pull/299)).
- Preserved per-hit truth during event overlays and derived optical-model bar lengths and readout positions from each bar's geometry ([#301](https://github.com/DUNE/dune-tms/pull/301)).
- Initialized the Kalman-node chi-squared values, recorded a covariance-weighted chi-squared, and corrected the muon/antimuon fit association ([#305](https://github.com/DUNE/dune-tms/pull/305), [#309](https://github.com/DUNE/dune-tms/pull/309)).
- Prevented index underflow and infinite loops while matching X/Y and U/V(/X) track candidates ([#310](https://github.com/DUNE/dune-tms/pull/310)).
- Corrected photon-timing simulation to sample individual scintillator and WLS delays and use the earliest detected-photon arrival, rather than an incorrectly parameterized Gamma shortcut ([#318](https://github.com/DUNE/dune-tms/pull/318)).

## [0.2] - 2023-12-13

- Areal density calculation fixed, affected Muon KE reconstruction
- Fixes to make larsoft work. Specifically scaling/units.
- Bugfix for infinite loop when a hit is outside volDetEnclosure
- Implemented simple timing simulation
- Added time slicer
- Added basic event display
- Added make_hists.py and draw_spill.py plotting scripts to provide examples in using the output
- Allow addressing bars as x, y, u, v; u/v are the +-3 degree planes
- Fix to optical fibre length in simulation
- Various minor bugfixes

## [TMSonlyFreeze] - 2022-03-01

- Tag for frozen TMS reconstruction. Used for the flux studies on TMS only

## [TrackMatching] - 2021-06-28

- Used for track matching studies with third production
