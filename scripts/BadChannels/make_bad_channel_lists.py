#!/usr/bin/env python3
"""Draw random bad channel lists for the bad channel scan (see src/TMS_BadChannels.h).

Input: a channel list from `DumpChannels` (CSV: plane,bar,view). For each seed the channels are put
in one fixed random order and a fraction f of them is the first round(f * N) of that order, so the
lists of one seed are nested (5% inside 10% inside ...) and the scan over f differs only by the
added channels, not by a new random draw. Different seeds are independent draws.

Output: <outdir>/f<percent:02d>_s<seed>.toml for each fraction and seed, plus f00.toml (no bad channels).
All channels are written with Status = "dead". Reproducible: same channels file, fractions and seeds
give identical files.

Usage: make_bad_channel_lists.py channels.csv outdir [--fractions 5 10 15 20 25] [--seeds 1 2 3]
"""
import argparse
import csv
import os
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("channels")
    ap.add_argument("outdir")
    ap.add_argument("--fractions", type=float, nargs="+", default=[5, 10, 15, 20, 25], help="percent of channels")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    args = ap.parse_args()

    with open(args.channels) as f:
        rows = list(csv.DictReader(f))
    channels = sorted((int(r["plane"]), int(r["bar"]), int(r["view"])) for r in rows)
    if len(set(channels)) != len(channels):
        raise SystemExit("channel list has duplicates")
    n = len(channels)
    os.makedirs(args.outdir, exist_ok=True)

    def write(name, chosen, note):
        with open(os.path.join(args.outdir, name), "w") as out:
            out.write(f"# {note}\n# Channels: {len(chosen)} of {n} ({100.0 * len(chosen) / n:.2f}%). See src/TMS_BadChannels.h.\n")
            out.write("BadChannel = [\n")
            for p, b, v in sorted(chosen):
                out.write(f'  {{Plane = {p}, Bar = {b}, View = {v}, Status = "dead"}},\n')
            out.write("]\n")

    write("f00.toml", [], "No bad channels.")
    for seed in args.seeds:
        order = list(channels)
        random.Random(seed).shuffle(order)
        for pct in args.fractions:
            k = round(pct / 100.0 * n)
            name = f"f{int(round(pct)):02d}_s{seed}.toml"
            write(name, order[:k], f"Random bad channels, {pct:g}% of channels, seed {seed} (nested within the larger fractions of this seed).")
            print(name, k)


if __name__ == "__main__":
    main()
