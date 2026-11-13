#!/usr/bin/env python3
"""Generate reproducible TMS field plots using the existing plot_mapper.py."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def grid_center(path):
    """Read grid origin, spacing and dimensions from mapper metadata."""
    import math

    origin = spacing = counts = None

    with path.open() as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                tokens = line[1:].strip().split()
                if tokens and tokens[0] == "shape":
                    tokens = tokens[1:]
                if len(tokens) >= 3:
                    try:
                        candidate = [int(v) for v in tokens[:3]]
                        if all(v > 0 for v in candidate):
                            counts = candidate
                    except ValueError:
                        pass
                continue

            values = line.split()
            if len(values) == 6:
                origin = [float(v) for v in values[:3]]
                spacing = [float(v) for v in values[3:]]
                break

    if origin is None or spacing is None:
        raise ValueError("Could not find mapper grid header")

    if counts is None:
        raise ValueError(
            "Grid dimensions not found in metadata; inspect header before plotting"
        )

    return {
        axis: origin[i] + spacing[i] * (counts[i] - 1) / 2
        for i, axis in enumerate("xyz")
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mapper",
        type=Path,
        default=ROOT / "Mapper_global_DIAGNOSTIC.txt"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "plots" / "automated_diagnostics"
    )
    parser.add_argument("--dpi", type=int, default=160)
    args = parser.parse_args()

    mapper = args.mapper.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    center = grid_center(mapper)
    commands = []

    for plane, orthogonal_axis in [
        ("xy", "z"),
        ("xz", "y"),
        ("yz", "x"),
    ]:
        for kind in ("2d", "quiver"):
            destination = output / f"field_{plane}_{kind}.png"

            command = [
                sys.executable,
                str(ROOT / "plot_mapper.py"),
                str(mapper),
                "--kind", kind,
                "--geometry", "stored",
                "--component", "mag",
                "--plane", plane,
                "--slice", str(center[orthogonal_axis]),
                "--units", "mm",
                "--dpi", str(args.dpi),
                "--output", str(destination),
            ]

            if kind == "quiver":
                command += ["--stride", "5"]

            print(f"Generating {destination.name}", flush=True)
            subprocess.run(command, check=True, cwd=ROOT)
            commands.append(command)

    manifest = {
        "mapper": str(mapper),
        "slice_centers_mm": center,
        "commands": commands,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(f"\nCompleted: {len(commands)} plots")
    print(f"Saved to: {output}")

if __name__ == "__main__":
    main()
