#!/usr/bin/env python3
"""
TMS_FLAT_CENTRAL_STEEL_REPRESENTATIVE_B.py

Canonical MAIN representative-field analysis for the DUNE TMS
FLAT CENTRAL STEEL REGION.

ANALYSIS CHAIN
--------------
This is main result #3 in the representative-field sequence:

  1. Full steel plate
  2. Central steel region
  3. Flat central steel region      <-- this script

Everything from the canonical central-steel analysis is kept fixed:
  - same 80 physical plates
  - same R1/R2/R3 plate-thickness geometry
  - same transverse central-steel selection
  - same physical thickness weighting
  - no |B| threshold
  - no smoothing

Only the longitudinal selection changes: within each plate-thickness family,
we retain the fixed z ranges chosen to represent the flatter/stable portions
of the central-steel field profile.

CENTRAL TRANSVERSE REGION
-------------------------
    10 <= |x| <= 3730 mm
    |y| <= 1780 mm

The 20 mm centre gap (-10 < x < +10 mm) is excluded.

PHYSICAL PLATE-THICKNESS REGIONS
--------------------------------
    R1 / Thin         : 34 plates x 15 mm
    R2 / Thick        : 22 plates x 40 mm
    R3 / Double-thick : 24 plates x 80 mm

FIXED FLAT LONGITUDINAL WINDOWS
-------------------------------
The current canonical local analysis uses plate-centre selections:

    R1: -3300 < z_center < -2300 mm
    R2: -1800 < z_center <   200 mm
    R3:   300 < z_center <  3000 mm

These ranges are intentionally explicit and are NOT found automatically.

METHOD
------
For every physical plate:
    B_i = mean |B| over central-steel mapper cells belonging to that plate

Select only plates whose physical centre lies in the fixed flat window
for their region.

For each selected region:
    B_rep(region) = sum_i(B_i * t_i) / sum_i(t_i)

For the combined flat-central selection:
    B_rep(flat central) = sum_i(B_i * t_i) / sum_i(t_i)

Because the same physical-thickness formula is used here as in the full-steel
and central-steel analyses, the only methodological change is the fixed
longitudinal selection.

Expected historical reference from the previous local analysis:
    R1 ~ 1.374 T
    R2 ~ 1.313 T
    R3 ~ 1.293 T
    Flat-central combined ~ 1.307 T

The script recomputes all values from Mapper.txt; none are hard-coded.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SHAPE_RE = re.compile(r"\bshape\s+(\d+)\s+(\d+)\s+(\d+)\b", re.I)

# ---------------------------------------------------------------------------
# Physical plate geometry from the current field dataset
# ---------------------------------------------------------------------------

THIN_STARTS_MM = np.arange(34, dtype=float) * 65.0
THICK_STARTS_MM = 2210.0 + np.arange(22, dtype=float) * 90.0
DOUBLE_STARTS_MM = 4190.0 + np.arange(24, dtype=float) * 135.0

SECTION_ORDER = ("R1 Thin", "R2 Thick", "R3 Double-thick")

SECTION_MARKERS = {
    "R1 Thin": "o",
    "R2 Thick": "s",
    "R3 Double-thick": "^",
}

# Fixed flat-region windows in TMS local z [mm].
FLAT_Z_WINDOWS = {
    "R1 Thin": (-3300.0, -2300.0),
    "R2 Thick": (-1800.0, 200.0),
    "R3 Double-thick": (300.0, 3000.0),
}


@dataclass(frozen=True)
class Plate:
    number: int
    section: str
    thickness_mm: float
    source_start_mm: float
    local_start_mm: float
    local_end_mm: float
    local_center_mm: float


# ---------------------------------------------------------------------------
# Mapper / geometry helpers
# ---------------------------------------------------------------------------

def read_mapper_header(path: Path):
    shape = None
    header = None

    with path.open("r") as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith("#"):
                m = SHAPE_RE.search(s)
                if m:
                    shape = tuple(int(v) for v in m.groups())
                continue

            values = s.split()

            if len(values) != 6:
                raise ValueError(
                    "Expected Mapper.txt header: ox oy oz dx dy dz"
                )

            header = tuple(float(v) for v in values)
            break

    if header is None:
        raise ValueError("Could not read Mapper.txt grid header.")

    return shape, header


def make_plates(local_z_origin_mm: float):
    plates = []
    number = 1

    definitions = (
        ("R1 Thin", 15.0, THIN_STARTS_MM),
        ("R2 Thick", 40.0, THICK_STARTS_MM),
        ("R3 Double-thick", 80.0, DOUBLE_STARTS_MM),
    )

    for section, thickness, starts in definitions:
        for source_start in starts:
            z0 = local_z_origin_mm + float(source_start)
            z1 = z0 + thickness

            plates.append(
                Plate(
                    number=number,
                    section=section,
                    thickness_mm=thickness,
                    source_start_mm=float(source_start),
                    local_start_mm=z0,
                    local_end_mm=z1,
                    local_center_mm=0.5 * (z0 + z1),
                )
            )

            number += 1

    if len(plates) != 80:
        raise RuntimeError(
            f"Expected 80 physical plates, found {len(plates)}."
        )

    expected_counts = {
        "R1 Thin": 34,
        "R2 Thick": 22,
        "R3 Double-thick": 24,
    }

    counts = {
        section: sum(p.section == section for p in plates)
        for section in SECTION_ORDER
    }

    if counts != expected_counts:
        raise RuntimeError(
            f"Unexpected plate-family counts: {counts}"
        )

    total_steel_mm = sum(
        p.thickness_mm
        for p in plates
    )

    if not math.isclose(total_steel_mm, 3310.0):
        raise RuntimeError(
            f"Expected 3310 mm total steel thickness, got {total_steel_mm} mm."
        )

    return plates


def full_steel_xy_mask(x: float, y: float) -> bool:
    """
    Same first-pass full-steel mask used in the canonical full-steel and
    central-steel analyses.
    """
    qx = abs(x)
    qy = abs(y)

    if qx < 10.0 or qx > 3730.0:
        return False

    if qy > 2350.0:
        return False

    if qy > 1780.0:
        if qx < 250.0 or qx > 3480.0:
            return False

    if 1620.0 <= qx <= 2120.0 and 1780.0 <= qy <= 2000.0:
        return False

    return True


def central_steel_xy_mask(x: float, y: float) -> bool:
    """
    Canonical central-steel transverse selection.

        10 <= |x| <= 3730 mm
        |y| <= 1780 mm
    """
    if not full_steel_xy_mask(x, y):
        return False

    qx = abs(x)
    qy = abs(y)

    return (
        10.0 <= qx <= 3730.0
        and qy <= 1780.0
    )


def build_z_to_plate(plates, z_values, tol=1e-6):
    """
    Assign each regular mapper z plane to a physical steel plate.

    Planes in longitudinal air gaps are mapped to None.
    """
    mapping = {}

    for z in z_values:
        idx = None

        for i, plate in enumerate(plates):
            if (
                z >= plate.local_start_mm - tol
                and z <= plate.local_end_mm + tol
            ):
                idx = i
                break

        mapping[round(float(z), 6)] = idx

    return mapping


def plate_is_in_flat_window(plate: Plate) -> bool:
    lo, hi = FLAT_Z_WINDOWS[plate.section]

    # Strict inequalities intentionally preserve the previous canonical script.
    return (
        lo < plate.local_center_mm < hi
    )


# ---------------------------------------------------------------------------
# Streaming accumulation over central steel
# ---------------------------------------------------------------------------

def new_accumulator(n):
    return {
        "count": np.zeros(n, dtype=np.int64),
        "sum_bmag": np.zeros(n, dtype=np.float64),
    }


def accumulate(acc, idx, bmag):
    acc["count"][idx] += 1
    acc["sum_bmag"][idx] += bmag


def stream_mapper(mapper_path, shape, header, plates):
    _, _, oz, _, _, dz = header

    if shape is None:
        max_end = max(
            p.local_end_mm
            for p in plates
        )

        nz = int(
            math.ceil(
                (max_end - oz) / dz
            )
        ) + 1
    else:
        _, _, nz = shape

    z_values = (
        oz
        + np.arange(
            nz,
            dtype=float,
        ) * dz
    )

    z_to_plate = build_z_to_plate(
        plates,
        z_values,
    )

    central_acc = new_accumulator(
        len(plates)
    )

    header_seen = False
    data_rows = 0
    plate_rows = 0
    central_rows = 0

    with mapper_path.open("r") as f:
        for line in f:
            s = line.strip()

            if not s or s.startswith("#"):
                continue

            values = s.split()

            if not header_seen:
                if len(values) != 6:
                    raise ValueError(
                        "Expected six-number Mapper.txt header."
                    )

                header_seen = True
                continue

            if len(values) < 6:
                continue

            x, y, z, bx, by, bz = map(
                float,
                values[:6],
            )

            if len(values) >= 7:
                bmag = float(
                    values[6]
                )
            else:
                bmag = math.sqrt(
                    bx * bx
                    + by * by
                    + bz * bz
                )

            data_rows += 1

            plate_idx = z_to_plate.get(
                round(z, 6)
            )

            if plate_idx is None:
                continue

            plate_rows += 1

            if central_steel_xy_mask(x, y):
                central_rows += 1

                accumulate(
                    central_acc,
                    plate_idx,
                    bmag,
                )

    return (
        central_acc,
        data_rows,
        plate_rows,
        central_rows,
    )


# ---------------------------------------------------------------------------
# Plate / region summaries
# ---------------------------------------------------------------------------

def build_plate_rows(plates, acc):
    rows = []

    for i, plate in enumerate(plates):
        n = int(
            acc["count"][i]
        )

        mean_b = (
            float(
                acc["sum_bmag"][i] / n
            )
            if n > 0
            else float("nan")
        )

        rows.append(
            {
                "plate": plate.number,
                "section": plate.section,
                "thickness_mm": plate.thickness_mm,
                "source_start_mm": plate.source_start_mm,
                "local_start_mm": plate.local_start_mm,
                "local_end_mm": plate.local_end_mm,
                "local_center_mm": plate.local_center_mm,
                "n_cells": n,
                "mean_B_T": mean_b,
                "selected_flat_window": int(
                    plate_is_in_flat_window(
                        plate
                    )
                ),
            }
        )

    return rows


def thickness_weighted_mean(rows):
    selected = [
        r
        for r in rows
        if np.isfinite(
            r["mean_B_T"]
        )
    ]

    if not selected:
        raise RuntimeError(
            "No valid rows available for thickness-weighted mean."
        )

    denominator = sum(
        r["thickness_mm"]
        for r in selected
    )

    numerator = sum(
        r["mean_B_T"]
        * r["thickness_mm"]
        for r in selected
    )

    return numerator / denominator


def summarize_region(selected_rows, section):
    rows = [
        r
        for r in selected_rows
        if r["section"] == section
    ]

    if not rows:
        raise RuntimeError(
            f"No selected flat-region plates found for {section}."
        )

    lo, hi = FLAT_Z_WINDOWS[section]

    return {
        "section": section,
        "z_min_mm": lo,
        "z_max_mm": hi,
        "n_plates": len(rows),
        "plate_thickness_mm": rows[0]["thickness_mm"],
        "selected_steel_thickness_mm": sum(
            r["thickness_mm"]
            for r in rows
        ),
        "Brep_T": thickness_weighted_mean(
            rows
        ),
    }


def summarize_flat_selection(plate_rows):
    selected = [
        r
        for r in plate_rows
        if (
            r["selected_flat_window"] == 1
            and np.isfinite(
                r["mean_B_T"]
            )
        )
    ]

    sections = [
        summarize_region(
            selected,
            section,
        )
        for section in SECTION_ORDER
    ]

    full = {
        "section": "Flat central steel",
        "n_plates": len(selected),
        "selected_steel_thickness_mm": sum(
            r["thickness_mm"]
            for r in selected
        ),
        "Brep_T": thickness_weighted_mean(
            selected
        ),
    }

    return selected, sections, full


def summarize_all_central(plate_rows):
    selected = [
        r
        for r in plate_rows
        if np.isfinite(
            r["mean_B_T"]
        )
    ]

    return {
        "section": "All central steel",
        "n_plates": len(selected),
        "steel_thickness_mm": sum(
            r["thickness_mm"]
            for r in selected
        ),
        "Brep_T": thickness_weighted_mean(
            selected
        ),
    }


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(path, rows):
    if not rows:
        return

    with path.open(
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(
                rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(rows)


def plot_flat_central_profile(
    path,
    plate_rows,
    sections,
    flat_summary,
    central_summary,
):
    fig, ax = plt.subplots(
        figsize=(13.5, 7.0)
    )

    # Background: all 80 central-steel physical plate means.
    ax.plot(
        [
            r["local_center_mm"]
            for r in plate_rows
        ],
        [
            r["mean_B_T"]
            for r in plate_rows
        ],
        marker="o",
        markersize=4,
        linewidth=1.0,
        alpha=0.25,
        label=(
            "All central-steel physical plate means "
            f"(B_rep={central_summary['Brep_T']:.3f} T)"
        ),
    )

    # Highlight fixed flat selections region by region.
    for sec in sections:
        section = sec["section"]

        rows = [
            r
            for r in plate_rows
            if (
                r["section"] == section
                and r["selected_flat_window"] == 1
            )
        ]

        z = [
            r["local_center_mm"]
            for r in rows
        ]

        b = [
            r["mean_B_T"]
            for r in rows
        ]

        ax.plot(
            z,
            b,
            marker=SECTION_MARKERS[section],
            markersize=5,
            linewidth=1.8,
            label=(
                f"{section} flat window: "
                f"B_rep={sec['Brep_T']:.3f} T"
            ),
        )

        ax.hlines(
            sec["Brep_T"],
            sec["z_min_mm"],
            sec["z_max_mm"],
            linestyles="--",
            linewidth=1.4,
        )

        ax.axvspan(
            sec["z_min_mm"],
            sec["z_max_mm"],
            alpha=0.06,
        )

    ax.axhline(
        flat_summary["Brep_T"],
        linestyle="-.",
        linewidth=2.0,
        label=(
            "Flat-central thickness-weighted "
            f"B_rep={flat_summary['Brep_T']:.3f} T"
        ),
    )

    ax.set_xlabel(
        "TMS local Z [mm]"
    )

    ax.set_ylabel(
        "Per-plate central-steel mean |B| [T]"
    )

    ax.set_title(
        "TMS Representative Magnetic Field — Flat Central Steel Region\n"
        "Same central transverse selection; fixed longitudinal flat windows"
    )

    ax.text(
        0.01,
        0.02,
        (
            "Fixed windows: "
            "R1 (-3300,-2300), "
            "R2 (-1800,200), "
            "R3 (300,3000) mm"
        ),
        transform=ax.transAxes,
        fontsize=9,
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        fontsize=8.5,
        loc="best",
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_central_vs_flat(
    path,
    central_summary,
    sections,
    flat_summary,
):
    labels = [
        "Central steel\n(all plates)",
        "Flat central\n(selected windows)",
    ]

    values = [
        central_summary["Brep_T"],
        flat_summary["Brep_T"],
    ]

    fig, ax = plt.subplots(
        figsize=(7.8, 6.0)
    )

    bars = ax.bar(
        labels,
        values,
    )

    for bar, value in zip(
        bars,
        values,
    ):
        ax.text(
            bar.get_x()
            + bar.get_width() / 2,
            value,
            f"{value:.4f} T",
            ha="center",
            va="bottom",
            fontsize=11,
        )

    ax.set_ylabel(
        "Thickness-weighted representative |B| [T]"
    )

    ax.set_title(
        "Central Steel vs Flat Central Steel\n"
        "Only the longitudinal selection changes"
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def write_summary(
    path,
    sections,
    flat_summary,
    central_summary,
    data_rows,
    plate_rows_count,
    central_rows,
):
    lines = [
        "TMS FLAT-CENTRAL-STEEL REPRESENTATIVE-FIELD ANALYSIS",
        "=" * 56,
        "",
        "MAIN ANALYSIS CHAIN",
        "  Full steel plate -> Central steel region -> Flat central steel region",
        "",
        "THIS SCRIPT",
        "  Main result #3: flat central steel region",
        "",
        "CENTRAL TRANSVERSE REGION",
        "  10 <= |x| <= 3730 mm",
        "  |y| <= 1780 mm",
        "  20 mm centre gap excluded",
        "",
        "PHYSICAL PLATE-THICKNESS REGIONS",
        "  R1 Thin         : 34 x 15 mm",
        "  R2 Thick        : 22 x 40 mm",
        "  R3 Double-thick : 24 x 80 mm",
        "",
        "FIXED FLAT LONGITUDINAL WINDOWS",
        "  R1: -3300 < z_center < -2300 mm",
        "  R2: -1800 < z_center <   200 mm",
        "  R3:   300 < z_center <  3000 mm",
        "",
        "SELECTION",
        "  Geometry only",
        "  No |B| threshold",
        "  No smoothing",
        "  Flat windows are fixed explicitly, not found automatically",
        "",
        "METHOD",
        "  For each physical plate:",
        "      B_i = mean |B| over central-steel cells",
        "",
        "  Keep only plates whose physical centre lies in the fixed",
        "  flat window for its thickness region.",
        "",
        "  B_rep = sum_i(B_i * t_i) / sum_i(t_i)",
        "",
        "FLAT-CENTRAL RESULTS",
    ]

    for sec in sections:
        lines.append(
            f"  {sec['section']:<18} "
            f"z=({sec['z_min_mm']:.0f},{sec['z_max_mm']:.0f}) mm  "
            f"n={sec['n_plates']:2d}  "
            f"selected steel={sec['selected_steel_thickness_mm']:.0f} mm  "
            f"B_rep={sec['Brep_T']:.6f} T"
        )

    lines.extend(
        [
            "",
            f"  Flat-central combined B_rep={flat_summary['Brep_T']:.6f} T",
            f"  Selected plates={flat_summary['n_plates']}",
            (
                "  Selected steel thickness="
                f"{flat_summary['selected_steel_thickness_mm']:.0f} mm"
            ),
            "",
            "CENTRAL -> FLAT-CENTRAL COMPARISON",
            f"  Central steel (all 80 plates): {central_summary['Brep_T']:.6f} T",
            f"  Flat central selection       : {flat_summary['Brep_T']:.6f} T",
            "",
            "ROW COUNTS",
            f"  Mapper data rows: {data_rows:,}",
            f"  Rows inside physical plate z spans: {plate_rows_count:,}",
            f"  Central-steel rows: {central_rows:,}",
            "",
            "IMPORTANT CAVEATS",
            "  * These flat z windows are analysis choices and should remain",
            "    explicit in plots, summaries, slides, and the repository.",
            "  * The 80 mm source-map offsets use 135 mm pitch; the current",
            "    simulation GDML has been observed to use 130 mm.",
            "  * This is a main representative-field result, distinct from the",
            "    automatic physical-plate cross-check and Maxwell/AEDT validation.",
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Canonical TMS flat-central-steel representative-field analysis."
        )
    )

    parser.add_argument(
        "mapper",
        nargs="?",
        default="Mapper.txt",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="plots/TMS_flat_central_steel_representative_B",
    )

    args = parser.parse_args()

    mapper_path = Path(
        args.mapper
    )

    outdir = Path(
        args.output_dir
    )

    if not mapper_path.exists():
        raise FileNotFoundError(
            mapper_path
        )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    shape, header = read_mapper_header(
        mapper_path
    )

    _, _, oz, _, _, _ = header

    plates = make_plates(
        oz
    )

    print(
        "TMS FLAT-CENTRAL-STEEL REPRESENTATIVE-FIELD ANALYSIS"
    )

    print(
        "====================================================="
    )

    print()

    print(
        "Central transverse selection:"
    )

    print(
        "  10 <= |x| <= 3730 mm"
    )

    print(
        "  |y| <= 1780 mm"
    )

    print(
        "  20 mm centre gap excluded"
    )

    print()

    print(
        "Fixed longitudinal flat windows:"
    )

    print(
        "  R1: -3300 < z_center < -2300 mm"
    )

    print(
        "  R2: -1800 < z_center <   200 mm"
    )

    print(
        "  R3:   300 < z_center <  3000 mm"
    )

    print()

    print(
        "Primary definition:"
    )

    print(
        "  same physical plates, transverse selection, and thickness weighting"
    )

    print(
        "  as the central-steel analysis"
    )

    print(
        "  only the longitudinal plate selection changes"
    )

    print(
        "  no B threshold"
    )

    print(
        "  no smoothing"
    )

    print()

    print(
        "Reading Mapper.txt in one pass..."
    )

    (
        central_acc,
        data_rows,
        plate_rows_count,
        central_rows,
    ) = stream_mapper(
        mapper_path,
        shape,
        header,
        plates,
    )

    plate_results = build_plate_rows(
        plates,
        central_acc,
    )

    (
        flat_rows,
        flat_sections,
        flat_summary,
    ) = summarize_flat_selection(
        plate_results
    )

    central_summary = summarize_all_central(
        plate_results
    )

    # Save all 80 plates with a selected/not-selected flag.
    write_csv(
        outdir / "flat_central_all_plate_results.csv",
        plate_results,
    )

    # Save only selected flat plates.
    write_csv(
        outdir / "flat_central_selected_plate_results.csv",
        flat_rows,
    )

    # Save region summaries plus combined summary.
    write_csv(
        outdir / "flat_central_section_results.csv",
        flat_sections + [flat_summary],
    )

    plot_flat_central_profile(
        outdir / "01_flat_central_steel_plate_profile.png",
        plate_results,
        flat_sections,
        flat_summary,
        central_summary,
    )

    plot_central_vs_flat(
        outdir / "02_central_vs_flat_central_comparison.png",
        central_summary,
        flat_sections,
        flat_summary,
    )

    write_summary(
        outdir / "flat_central_steel_representative_B_summary.txt",
        flat_sections,
        flat_summary,
        central_summary,
        data_rows,
        plate_rows_count,
        central_rows,
    )

    print()

    print(
        "FLAT-CENTRAL MAIN RESULTS"
    )

    print(
        "-------------------------"
    )

    for sec in flat_sections:
        print(
            f"{sec['section']:<18} "
            f"z=({sec['z_min_mm']:.0f},{sec['z_max_mm']:.0f}) mm  "
            f"n={sec['n_plates']:2d}  "
            f"B_rep={sec['Brep_T']:.4f} T"
        )

    print()

    print(
        "Flat-central thickness-weighted "
        f"B_rep = {flat_summary['Brep_T']:.4f} T"
    )

    print()

    print(
        "CENTRAL STEEL -> FLAT CENTRAL STEEL"
    )

    print(
        "-----------------------------------"
    )

    change = (
        100.0
        * (
            central_summary["Brep_T"]
            - flat_summary["Brep_T"]
        )
        / central_summary["Brep_T"]
    )

    print(
        f"Central steel     : {central_summary['Brep_T']:.6f} T"
    )

    print(
        f"Flat central steel: {flat_summary['Brep_T']:.6f} T"
    )

    print(
        f"Change            : {change:.2f}% lower"
    )

    print()

    print(
        "Saved:"
    )

    print(
        f"  {outdir / '01_flat_central_steel_plate_profile.png'}"
    )

    print(
        f"  {outdir / '02_central_vs_flat_central_comparison.png'}"
    )

    print(
        f"  {outdir / 'flat_central_all_plate_results.csv'}"
    )

    print(
        f"  {outdir / 'flat_central_selected_plate_results.csv'}"
    )

    print(
        f"  {outdir / 'flat_central_section_results.csv'}"
    )

    print(
        f"  {outdir / 'flat_central_steel_representative_B_summary.txt'}"
    )


if __name__ == "__main__":
    main()
