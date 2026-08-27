# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to a modified version of [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Geometry releases will be tagged as `Descriptive_tag_v_X.Y.Z`.

## [1.1.0] - Unreleased

### Added

- Added run-aware global vertex identifiers for truth tracking ([#247](https://github.com/DUNE/dune-tms/pull/247)).
- Restored the GitHub Actions build workflow ([#259](https://github.com/DUNE/dune-tms/pull/259)).
- Added `nTrueParticles` back to `Truth_Info` ([#282](https://github.com/DUNE/dune-tms/pull/282)).
- Added the missing Y-view line-candidate output and corrected X/Y branch names ([#284](https://github.com/DUNE/dune-tms/pull/284)).
- Added leading truth-contributor information for reconstructed hits ([#303](https://github.com/DUNE/dune-tms/pull/303)).

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
- Aligned the CMake package version and ROOT output metadata version at `1.1.0`.
- Changed `MergeTracks` to default to true.
- Enabled adding trailing Hough-track hits to their corresponding track by default, with a configuration switch.

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
- Prevented 3D track matching from looping indefinitely when the last hit in a required view lies outside the detector bounds ([#310](https://github.com/DUNE/dune-tms/pull/310)).

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
