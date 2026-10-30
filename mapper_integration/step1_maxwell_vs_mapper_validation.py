#!/usr/bin/env python3
"""
Maxwell/AEDT vs uniform Mapper validation for the DUNE TMS field map.

Outputs:
  - four individual Maxwell-vs-Mapper longitudinal profile plots
  - matched-point CSVs for P1-P4
  - validation_metrics.csv
  - one combined 2x2 figure with all four polylines

Notes:
  - Maxwell/AEDT validation uses |B| > 0.5 T.
  - 10-point rolling mean is for DISPLAY only.
  - Numerical metrics are calculated from unsmoothed matched points.
  - P1/P4 use Mapper x=100 mm as the nearest steel-side proxy
    for Maxwell x=10 mm.
  - P2/P3 are nearly coordinate matched:
    Maxwell x=2795 mm vs Mapper x=2800 mm.
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


POLYLINES = {
    1: {"file": "polyline1_30_400.csv", "x": 10.0, "y": 0.0},
    2: {"file": "polyline2_30_400.csv", "x": 2795.0, "y": 0.0},
    3: {"file": "polyline3_30_400.csv", "x": 2795.0, "y": 900.0},
    4: {"file": "polyline4_30_400.csv", "x": 10.0, "y": 900.0},
}

Z_MIN = -4000.0
Z_MAX = 3260.0
B_THRESHOLD = 0.5
ROLLING_WINDOW = 10


def read_mapper_header(mapper_path):
    with open(mapper_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            vals = s.split()
            if len(vals) != 6:
                raise ValueError("Expected Mapper header: ox oy oz dx dy dz")
            return tuple(map(float, vals))
    raise RuntimeError("Could not read Mapper.txt header.")


def choose_mapper_xy(ox, oy, dx, dy, x_req, y_req):
    if abs(x_req - 10.0) < 1e-9:
        n = math.ceil((10.0 - ox) / dx)
        x_map = ox + n * dx
    else:
        x_map = ox + round((x_req - ox) / dx) * dx

    y_map = oy + round((y_req - oy) / dy) * dy
    return float(x_map), float(y_map)


def read_maxwell_csv(csv_path):
    df = pd.read_csv(csv_path)

    if "Distance [meter]" not in df.columns:
        raise ValueError(f"{csv_path.name}: missing 'Distance [meter]'")
    if "Mag_B [tesla]" not in df.columns:
        raise ValueError(f"{csv_path.name}: missing 'Mag_B [tesla]'")

    out = pd.DataFrame()
    out["distance_mm"] = df["Distance [meter]"].astype(float) * 1000.0
    out["z_mm"] = Z_MIN + out["distance_mm"]
    out["B_T"] = df["Mag_B [tesla]"].astype(float)

    return out[
        (out["z_mm"] >= Z_MIN)
        & (out["z_mm"] <= Z_MAX)
    ].copy()


def extract_mapper_profile(mapper_path, x_target, y_target):
    rows = []
    header_seen = False

    with open(mapper_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            vals = s.split()

            if not header_seen:
                header_seen = True
                continue

            if len(vals) < 6:
                continue

            x, y, z, bx, by, bz = map(float, vals[:6])

            if abs(x - x_target) > 1e-6:
                continue
            if abs(y - y_target) > 1e-6:
                continue
            if z < Z_MIN or z > Z_MAX:
                continue

            if len(vals) >= 7:
                bmag = float(vals[6])
            else:
                bmag = math.sqrt(bx * bx + by * by + bz * bz)

            rows.append(
                {
                    "z_mm": z,
                    "distance_mm": z - Z_MIN,
                    "B_T": bmag,
                }
            )

    if not rows:
        raise RuntimeError(
            f"No Mapper data found at x={x_target}, y={y_target}"
        )

    return (
        pd.DataFrame(rows)
        .sort_values("z_mm")
        .reset_index(drop=True)
    )


def process_profile(df):
    selected = df[df["B_T"] > B_THRESHOLD].copy()

    selected["B_smooth_T"] = (
        selected["B_T"]
        .rolling(
            ROLLING_WINDOW,
            center=True,
        )
        .mean()
    )

    return selected


def compare_profiles(maxwell_df, mapper_df):
    max_df = (
        maxwell_df[["z_mm", "B_T"]]
        .rename(columns={"B_T": "B_Maxwell_T"})
        .sort_values("z_mm")
    )

    map_df = (
        mapper_df[["z_mm", "B_T"]]
        .rename(columns={"B_T": "B_Mapper_T"})
        .sort_values("z_mm")
    )

    merged = pd.merge_asof(
        max_df,
        map_df,
        on="z_mm",
        direction="nearest",
        tolerance=6.0,
    ).dropna().copy()

    if len(merged) == 0:
        return merged, {}

    diff = merged["B_Mapper_T"] - merged["B_Maxwell_T"]

    metrics = {
        "N matched": len(merged),
        "Maxwell mean [T]": merged["B_Maxwell_T"].mean(),
        "Mapper mean [T]": merged["B_Mapper_T"].mean(),
        "Bias Mapper-Maxwell [T]": diff.mean(),
        "MAE [T]": diff.abs().mean(),
        "RMSE [T]": np.sqrt(np.mean(diff**2)),
        "Correlation": merged[
            ["B_Maxwell_T", "B_Mapper_T"]
        ].corr().iloc[0, 1],
    }

    return merged, metrics


def make_plot(
    polyline_number,
    maxwell_df,
    mapper_df,
    maxwell_xy,
    mapper_xy,
    outpath,
):
    fig, ax = plt.subplots(figsize=(11, 6))

    ax.plot(
        maxwell_df["distance_mm"],
        maxwell_df["B_smooth_T"],
        linewidth=2.0,
        label=f"Maxwell ({maxwell_xy[0]:g}, {maxwell_xy[1]:g}) mm",
    )

    ax.plot(
        mapper_df["distance_mm"],
        mapper_df["B_smooth_T"],
        linestyle="--",
        linewidth=2.0,
        label=f"Mapper ({mapper_xy[0]:g}, {mapper_xy[1]:g}) mm",
    )

    role = (
        "steel-side proxy"
        if polyline_number in (1, 4)
        else "nearly coordinate-matched"
    )

    ax.set_xlabel("Distance along TMS / beam direction [mm]")
    ax.set_ylabel("|B| [T]")
    ax.set_title(
        f"Polyline {polyline_number}: Maxwell/AEDT vs Mapper — {role}\n"
        f"Validation subset: |B| > {B_THRESHOLD} T; "
        f"{ROLLING_WINDOW}-point rolling mean shown only for visualization"
    )

    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_four_panel_plot(plot_payload, outpath):
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(16, 11),
        sharex=True,
        sharey=False,
    )
    axes = axes.flatten()

    for ax, number in zip(axes, sorted(plot_payload)):
        item = plot_payload[number]

        maxwell_df = item["maxwell"]
        mapper_df = item["mapper"]
        maxwell_xy = item["maxwell_xy"]
        mapper_xy = item["mapper_xy"]

        ax.plot(
            maxwell_df["distance_mm"],
            maxwell_df["B_smooth_T"],
            linewidth=2.0,
            label=f"Maxwell ({maxwell_xy[0]:g}, {maxwell_xy[1]:g}) mm",
        )

        ax.plot(
            mapper_df["distance_mm"],
            mapper_df["B_smooth_T"],
            linestyle="--",
            linewidth=2.0,
            label=f"Mapper ({mapper_xy[0]:g}, {mapper_xy[1]:g}) mm",
        )

        role = (
            "steel-side proxy"
            if number in (1, 4)
            else "nearly coordinate-matched"
        )

        ax.set_title(f"P{number} — {role}")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=9)

    axes[0].set_ylabel("|B| [T]")
    axes[2].set_ylabel("|B| [T]")
    axes[2].set_xlabel("Distance along TMS / beam direction [mm]")
    axes[3].set_xlabel("Distance along TMS / beam direction [mm]")

    fig.suptitle(
        "Maxwell/AEDT vs Mapper Validation — All Four Polylines\n"
        f"Validation subset: |B| > {B_THRESHOLD} T; "
        f"{ROLLING_WINDOW}-point rolling mean shown only for visualization",
        fontsize=16,
    )

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "mapper",
        help="Path to Mapper.txt",
    )

    parser.add_argument(
        "maxwell_dir",
        help="Folder containing Marco's four CSV files",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="plots/TMS_maxwell_vs_mapper_validation",
    )

    args = parser.parse_args()

    mapper_path = Path(args.mapper)
    maxwell_dir = Path(args.maxwell_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    ox, oy, oz, dx, dy, dz = read_mapper_header(mapper_path)

    print()
    print("MAXWELL/AEDT-vs-MAPPER VALIDATION")
    print("=================================")
    print()
    print(f"Mapper dx = {dx:g} mm")
    print(f"Mapper dy = {dy:g} mm")
    print(f"Mapper dz = {dz:g} mm")
    print()
    print(f"Validation subset: |B| > {B_THRESHOLD} T")
    print(f"Display smoothing: {ROLLING_WINDOW}-point rolling mean")
    print("Numerical metrics: unsmoothed matched points")
    print()

    all_metrics = []
    plot_payload = {}

    for number, info in POLYLINES.items():
        csv_path = maxwell_dir / info["file"]

        if not csv_path.exists():
            raise FileNotFoundError(csv_path)

        x_map, y_map = choose_mapper_xy(
            ox,
            oy,
            dx,
            dy,
            info["x"],
            info["y"],
        )

        role = (
            "steel-side proxy"
            if number in (1, 4)
            else "nearly coordinate-matched"
        )

        print(
            f"P{number}: Maxwell ({info['x']:g}, {info['y']:g}) mm "
            f"-> Mapper ({x_map:g}, {y_map:g}) mm [{role}]"
        )

        maxwell_raw = read_maxwell_csv(csv_path)

        mapper_raw = extract_mapper_profile(
            mapper_path,
            x_map,
            y_map,
        )

        maxwell = process_profile(maxwell_raw)
        mapper = process_profile(mapper_raw)

        merged, metrics = compare_profiles(
            maxwell,
            mapper,
        )

        if metrics:
            row = {
                "Polyline": number,
                "Maxwell x [mm]": info["x"],
                "Maxwell y [mm]": info["y"],
                "Mapper x [mm]": x_map,
                "Mapper y [mm]": y_map,
                "Role": role,
                **metrics,
            }
            row["Relative bias [%]"] = (
                100.0
                * row["Bias Mapper-Maxwell [T]"]
                / row["Maxwell mean [T]"]
            )
            all_metrics.append(row)

        merged.to_csv(
            outdir / f"polyline{number}_matched_points.csv",
            index=False,
        )

        make_plot(
            number,
            maxwell,
            mapper,
            (info["x"], info["y"]),
            (x_map, y_map),
            outdir / f"polyline{number}_maxwell_vs_mapper.png",
        )

        plot_payload[number] = {
            "maxwell": maxwell.copy(),
            "mapper": mapper.copy(),
            "maxwell_xy": (info["x"], info["y"]),
            "mapper_xy": (x_map, y_map),
        }

    make_four_panel_plot(
        plot_payload,
        outdir / "05_all_four_polylines_2x2.png",
    )

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(
        outdir / "validation_metrics.csv",
        index=False,
    )

    print()
    print("VALIDATION METRICS")
    print("==================")
    print()

    if len(metrics_df):
        print(
            metrics_df[
                [
                    "Polyline",
                    "Role",
                    "Maxwell mean [T]",
                    "Mapper mean [T]",
                    "Bias Mapper-Maxwell [T]",
                    "Relative bias [%]",
                    "MAE [T]",
                    "RMSE [T]",
                    "Correlation",
                ]
            ].to_string(index=False)
        )

    print()
    print("Important:")
    print(
        "P1/P4 use Mapper x=100 mm as the nearest steel-side proxy "
        "for Maxwell x=10 mm."
    )
    print(
        "P2/P3 use Mapper x=2800 mm for Maxwell x=2795 mm."
    )
    print()
    print(f"Results saved in: {outdir}")
    print(
        f"Combined 2x2 figure: "
        f"{outdir / '05_all_four_polylines_2x2.png'}"
    )


if __name__ == "__main__":
    main()
