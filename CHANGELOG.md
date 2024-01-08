# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to a modified version of [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Geometry releases will be tagged as `Descriptive_tag_v_X.Y.Z`.

## [Version 0.2] - 2023-12-13
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
