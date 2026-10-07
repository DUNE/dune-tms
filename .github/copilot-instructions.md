# Review instructions for dune-tms

dune-tms is the reconstruction and detector-simulation code for the DUNE Temporary Muon
Spectrometer (TMS). It is C++ built on ROOT, edep-sim and toml11, configured by TOML files.
Most of the value of a review here is catching bugs that change physics output silently, so
be specific and be sure.

## How to review

- **Read beyond the diff.** Many changes add a call to, or a use of, code in a file the PR
  does not touch. Before saying something is unused, never called, missing, or never reset,
  search the whole repository (`app/`, `src/`, `scripts/`, `config/`) and say what you searched.
- **Cite file and line for every finding**, and say what input or state makes it fail. Drop
  findings you cannot tie to a concrete failure.
- **Separate correctness from preference.** Report bugs first. Style, naming and refactoring
  suggestions go last and are labeled as optional. Do not request changes for them.
- **Do not repeat the PR description back** as a summary; say only what you found.
- If the PR says it was verified (for example "output identical, branch by branch"), review
  whether that check could have caught the problem, rather than assuming it did not run.

## What to check

- **Array and index bounds.** Output trees use fixed-size arrays (for example in
  `TMS_TreeWriter`); check that every fill is clamped or guarded. Check `-1` or empty-container
  indices before they are used, and signed/unsigned comparisons.
- **Units.** Positions are millimeters, times nanoseconds, energies MeV unless a name says
  otherwise. Flag mixed units, especially cm/mm and meters in Cluster3D.
- **Integer width.** Look for truncation of 64-bit values (event or spill numbers, seeds, IDs)
  into 32-bit types.
- **Determinism.** Output must not depend on pointer values, unordered-container iteration
  order, uninitialized memory, build flags, or the order files are processed. Random seeds
  must come from event identity, not wall-clock time.
- **Performance-only changes must not change output.** Reordering loops, caching, or
  replacing containers is fine only if results are bit-identical. Flag floating-point
  reassociation, changed iteration order, and changed tie-breaking.
- **Config keys.** New TOML keys should be read with a default (`find_or`-style), not a
  mandatory `toml::find` that breaks existing config files. Check the default in
  `config/TMS_Default_Config.toml` matches the one in code, and that a new key is documented
  there with a comment.
- **Output compatibility.** Renaming or removing a branch in `Reco_Tree`, `Truth_Info`,
  `Reco_Tree_C3D` or `Truth_Info_C3D` breaks the validation scripts in `scripts/Validation`.
  Flag it unless the PR updates the readers.
- **Memory and ownership.** Raw `new` without a matching delete, dangling references to
  vector elements after a push_back, and large copies in per-hit or per-event loops.
- **Geometry assumptions.** The detector has x- and y-measuring bars with different readout
  (vertical bars read at the top, horizontal bars split at x = 0 and read at the outer
  sides). Flag code that treats the two views as symmetric without saying so.

## Build and CI

- The project builds two ways: CMake (local and grid builds) and the plain `Makefile`s
  under `src/` and `app/`, which GitHub CI uses (`.github/workflows/c-cpp.yml`, pull
  requests to `main`). A new source file or include directory must be added to **both**, or
  CI fails with "no such file" or an unresolved symbol even when the CMake build passes.
- CMake defaults to a Debug build, and grid production uses an optimized one
  (`-DCMAKE_BUILD_TYPE=RelWithDebInfo`). Code must give the same results in both, so flag
  correctness that depends on `assert` or on uninitialized memory behaving a certain way.
- There is no unit-test suite. Changes are checked by running `ConvertToTMSTree` on
  production files and comparing output trees, and by the scripts in `scripts/Validation`.
  A PR that changes reconstruction behavior should say what it was compared against.

## Conventions

- American spelling (meter, center, millimeter) in comments, docs and messages.
- Comments should explain why, since the code is maintained by physicists as well as
  programmers. Prefer a short explanatory comment over none on non-obvious logic.
- Keep READMEs current: changes under `src/Cluster3D` should update `src/Cluster3D/README.md`.
- New behavior that changes physics output should be behind a config switch, with the
  default stated in the PR, so that results can be compared before and after.
- Python scripts should run on the Python 3.6 and newer found in the SL7 container.

## What not to flag

- `src/attic/` holds retired code that is not built.
- Do not ask for a unit-test framework, a formatter, or a rewrite of an existing file.
- Do not flag the unconditional `STAGETIME` summary printed by `ConvertToTMSTree`; it is
  intended output, used in grid logs.
- Third-party code in `toml11/` is out of scope.
