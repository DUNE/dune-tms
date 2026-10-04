#!/bin/bash
# Run the controlled single-stage scans behind the stage-by-stage detector-simulation page.
# Usage: run_stage_scans.sh <build dir> <edep-sim file for the geometry> <output dir> [throw scale]
# Needs the environment of setup.sh. Every run pins its configuration (copies of the repository defaults
# with the listed overrides, written to <output dir>/configs), since the defaults are otherwise read from
# $TMS_DIR at run time. A throw scale below 1 gives a quick test run.
set -e
B=$1; GEOM=$2; OUT=$3; SCALE=${4:-1}
[ -x "$B/app/DetSimStageScan" ] || { echo "no DetSimStageScan in $B/app"; exit 1; }
SRC=$(cd "$(dirname "$0")/../../../.." && pwd)
mkdir -p "$OUT/configs"
cp "$SRC/config/TMS_Default_Config.toml" "$OUT/configs/tms.toml"
n() { python3 -c "print(max(1,int($1*$SCALE)))"; }

# readout config = default + "Key = value" overrides (each key must exist exactly once in the default)
mkconf() {
  local name=$1; shift
  python3 - "$SRC/config/TMS_Readout_Default_Config.toml" "$OUT/configs/readout_$name.toml" "$@" <<'PY'
import re, sys
src, dst, *kv = sys.argv[1:]
text = open(src).read()
for item in kv:
    key, val = item.split("=", 1)
    pat = re.compile(r"^(\s*%s\s*=\s*)[^#\n]*?(\s*(#.*)?)$" % re.escape(key), re.M)
    assert len(pat.findall(text)) == 1, key
    text = pat.sub(lambda m: m.group(1) + val + m.group(2), text)
open(dst, "w").write(text)
PY
}
mkconf legacy UseResponseElements=false
mkconf new UseResponseElements=true
mkconf new_dead500 UseResponseElements=true Deadtime=500.0
mkconf new_dead500_zombie100 UseResponseElements=true Deadtime=500.0 ZombieTime=100.0
# deposit bin length 0.5 and 2 mm around the default 1 mm, for the response-grid convergence check
for b in 0.5 2.0; do mkconf new_bin$b UseResponseElements=true DepositBinLength=$b; done
THRS="0.5 1.0 1.5 2.0 2.5 3.0"
FAR=3430  # mm from the readout end of the 3.5 m reference bar: the far end
for t in $THRS; do mkconf timing_thr$t UseResponseElements=true FrontEndTimingMode=true DiscriminatorThreshold=$t; done

run() {  # config, log name, command...
  local conf=$1 name=$2; shift 2
  TMS_TOML=$OUT/configs/tms.toml TMS_READOUT_TOML=$OUT/configs/readout_$conf.toml "$@" > "$OUT/$name.log" 2>&1 \
    && echo "$name ok" || { echo "$name FAILED (see $OUT/$name.log)"; return 1; }
}
pids=()
for c in legacy new; do
  run $c path_far_$c $B/app/DetSimStageScan path "$GEOM" "$OUT/path_far_$c.csv" $(n 4000) $FAR & pids+=($!)
  run $c pileup_$c $B/app/DetSimStageScan pileup "$GEOM" "$OUT/pileup_$c.csv" $(n 400) & pids+=($!)
  run $c reseg_$c $B/app/ArtificialResegmentationTest "$GEOM" "$OUT/reseg_$c.csv" $(n 5000) & pids+=($!)
  run $c path_$c $B/app/DetSimStageScan path "$GEOM" "$OUT/path_$c.csv" $(n 4000) & pids+=($!)
  run $c position_$c $B/app/DetSimStageScan position "$GEOM" "$OUT/position_$c.csv" $(n 2000) & pids+=($!)
done
for c in legacy new new_dead500 new_dead500_zombie100; do
  run $c pair_$c $B/app/DetSimStageScan pair "$GEOM" "$OUT/pair_$c.csv" $(n 100) & pids+=($!)
done
for b in 0.5 2.0; do
  run new_bin$b path_new_bin$b $B/app/DetSimStageScan path "$GEOM" "$OUT/path_new_bin$b.csv" $(n 4000) & pids+=($!)
  run new_bin$b reseg_new_bin$b $B/app/ArtificialResegmentationTest "$GEOM" "$OUT/reseg_new_bin$b.csv" $(n 5000) & pids+=($!)
done
for t in $THRS; do
  run timing_thr$t path_timing_thr$t $B/app/DetSimStageScan path "$GEOM" "$OUT/path_timing_thr$t.csv" $(n 2000) & pids+=($!)
  run timing_thr$t path_timing_thr${t}_far $B/app/DetSimStageScan path "$GEOM" "$OUT/path_timing_thr${t}_far.csv" $(n 2000) $FAR & pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait $p || fail=1; done
( cd "$SRC" && echo "source: $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- src app config || echo ' + uncommitted changes')" ) > "$OUT/provenance.txt"
echo "build: $B" >> "$OUT/provenance.txt"; echo "geometry: $GEOM" >> "$OUT/provenance.txt"; echo "throw scale: $SCALE" >> "$OUT/provenance.txt"
exit $fail
