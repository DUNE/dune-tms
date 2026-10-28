#!/usr/bin/env python3
"""
TMS_CENTRAL_STEEL_REPRESENTATIVE_B.py

Canonical MAIN representative-field analysis for the DUNE TMS CENTRAL STEEL REGION.

ANALYSIS CHAIN
--------------
This is main result #2 in the representative-field sequence:

  1. Full steel plate
  2. Central steel region        <-- this script
  3. Flat central steel region

The physical longitudinal geometry is kept identical to the full-steel analysis:

    R1 / Thin         : 34 plates x 15 mm
    R2 / Thick        : 22 plates x 40 mm
    R3 / Double-thick : 24 plates x 80 mm

Only the transverse spatial selection is changed.

CENTRAL STEEL REGION
--------------------
The central transverse steel region is defined by

    10 <= |x| <= 3730 mm
    |y| <= 1780 mm

The 20 mm centre gap (-10 < x < +10 mm) is therefore excluded.

On the 100 mm mapper grid, the actual selected centres are approximately

    x = +/-100, +/-200, ... +/-3700 mm
    y = -1700, -1600, ... +1700 mm

No magnetic-field threshold is used.
No smoothing is used.

For each physical plate:

    B_i = mean |B| over selected central-steel mapper cells in that plate

For each thickness region:

    B_rep(region) = sum_i(B_i * t_i) / sum_i(t_i)

For the full central-steel detector selection:

    B_rep(central) = sum_i(B_i * t_i) / sum_i(t_i)

The script also recomputes the full-steel value with the same geometry-mask logic
for a direct full-steel -> central-steel comparison.

Expected historical reference for the current mapper:
    central full-TMS B_rep ~ 1.320 T

The numerical result is recomputed from Mapper.txt; it is not hard-coded.
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
# Mapper / plate geometry helpers
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

    counts = {
        section: sum(p.section == section for p in plates)
        for section in SECTION_ORDER
    }

    if counts != {
        "R1 Thin": 34,
        "R2 Thick": 22,
        "R3 Double-thick": 24,
    }:
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
    Same first-pass full-steel footprint used in the canonical full-steel analysis.
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

    The explicit full-steel-mask intersection documents that the central selection
    is a subset of the same steel geometry used by the full-steel analysis.
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


# ---------------------------------------------------------------------------
# Streaming accumulation
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

    full_acc = new_accumulator(
        len(plates)
    )

    central_acc = new_accumulator(
        len(plates)
    )

    header_seen = False
    data_rows = 0
    full_rows = 0
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

            if full_steel_xy_mask(x, y):
                full_rows += 1

                accumulate(
                    full_acc,
                    plate_idx,
                    bmag,
                )

            if central_steel_xy_mask(x, y):
                central_rows += 1

                accumulate(
                    central_acc,
                    plate_idx,
                    bmag,
                )

    return (
        full_acc,
        central_acc,
        data_rows,
        full_rows,
        central_rows,
    )


# ---------------------------------------------------------------------------
# Plate / section summaries
# ---------------------------------------------------------------------------

def plate_rows(plates, acc):
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
            }
        )

    return rows


def summarize_section(rows, section):
    selected = [
        r
        for r in rows
        if (
            r["section"] == section
            and np.isfinite(
                r["mean_B_T"]
            )
        )
    ]

    if not selected:
        raise RuntimeError(
            f"No valid plates found for {section}."
        )

    total_t = sum(
        r["thickness_mm"]
        for r in selected
    )

    brep = (
        sum(
            r["mean_B_T"]
            * r["thickness_mm"]
            for r in selected
        )
        / total_t
    )

    return {
        "section": section,
        "n_plates": len(selected),
        "plate_thickness_mm": selected[0]["thickness_mm"],
        "total_steel_thickness_mm": total_t,
        "Brep_T": brep,
    }


def summarize(rows):
    sections = [
        summarize_section(
            rows,
            section,
        )
        for section in SECTION_ORDER
    ]

    selected = [
        r
        for r in rows
        if np.isfinite(
            r["mean_B_T"]
        )
    ]

    total_t = sum(
        r["thickness_mm"]
        for r in selected
    )

    full = {
        "section": "Full TMS",
        "n_plates": len(selected),
        "plate_thickness_mm": "",
        "total_steel_thickness_mm": total_t,
        "Brep_T": (
            sum(
                r["mean_B_T"]
                * r["thickness_mm"]
                for r in selected
            )
            / total_t
        ),
    }

    return sections, full


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


def plot_region_definition(path):
    fig, ax = plt.subplots(
        figsize=(11.5, 6.3)
    )

    # Full simplified plate envelope.
    ax.add_patch(
        plt.Rectangle(
            (-3730, -2350),
            7460,
            4700,
            fill=False,
            linewidth=1.5,
            label="Full plate envelope",
        )
    )

    # Canonical central-steel selection.
    ax.add_patch(
        plt.Rectangle(
            (-3730, -1780),
            7460,
            3560,
            alpha=0.18,
            label="Selected central steel region",
        )
    )

    # 20 mm centre gap.
    ax.add_patch(
        plt.Rectangle(
            (-10, -1780),
            20,
            3560,
            facecolor="white",
            edgecolor="black",
            linewidth=1.0,
        )
    )

    ax.axhline(
        1780,
        linestyle=":",
        linewidth=1.1,
    )

    ax.axhline(
        -1780,
        linestyle=":",
        linewidth=1.1,
    )

    ax.axvline(
        0,
        linestyle="--",
        linewidth=1.0,
    )

    ax.text(
        0,
        0,
        "20 mm\ncentre gap",
        ha="center",
        va="center",
        fontsize=10,
    )

    ax.text(
        0,
        1500,
        "Selected central steel region",
        ha="center",
        va="center",
        fontsize=13,
    )

    ax.set_xlim(
        -4000,
        4000,
    )

    ax.set_ylim(
        -2600,
        2600,
    )

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.set_xlabel(
        "X [mm]"
    )

    ax.set_ylabel(
        "Y [mm]"
    )

    ax.set_title(
        "TMS Central Steel Region Used for Representative-Field Analysis\n"
        r"$10 \leq |x| \leq 3730$ mm,  $|y| \leq 1780$ mm"
    )

    ax.grid(
        True,
        alpha=0.15,
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_central_profile(path, rows, sections, full):
    fig, ax = plt.subplots(
        figsize=(13.2, 6.8)
    )

    for sec in sections:
        section = sec["section"]

        selected = [
            r
            for r in rows
            if r["section"] == section
        ]

        z = np.array(
            [
                r["local_center_mm"]
                for r in selected
            ]
        )

        b = np.array(
            [
                r["mean_B_T"]
                for r in selected
            ]
        )

        ax.plot(
            z,
            b,
            marker=SECTION_MARKERS[section],
            linewidth=1.4,
            markersize=5,
            label=(
                f"{section}: "
                f"B_rep={sec['Brep_T']:.3f} T"
            ),
        )

        ax.hlines(
            sec["Brep_T"],
            min(
                r["local_start_mm"]
                for r in selected
            ),
            max(
                r["local_end_mm"]
                for r in selected
            ),
            linestyles="--",
            linewidth=1.3,
        )

    ax.axhline(
        full["Brep_T"],
        linestyle="-.",
        linewidth=2.0,
        label=(
            "Central-steel thickness-weighted "
            f"B_rep={full['Brep_T']:.3f} T"
        ),
    )

    ax.set_xlabel(
        "TMS local Z [mm]"
    )

    ax.set_ylabel(
        "Per-plate mean |B| [T]"
    )

    ax.set_title(
        "TMS Representative Magnetic Field — Central Steel Region\n"
        r"$10 \leq |x| \leq 3730$ mm,  $|y| \leq 1780$ mm; "
        "same 15 / 40 / 80 mm plate-thickness regions"
    )

    ax.text(
        0.01,
        0.02,
        (
            "Primary central-steel result: geometry selection only, "
            "no B threshold, no smoothing"
        ),
        transform=ax.transAxes,
        fontsize=9,
    )

    ax.grid(
        True,
        alpha=0.22,
    )

    ax.legend(
        fontsize=9,
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_comparison(
    path,
    full_sections,
    central_sections,
    full_all,
    central_all,
):
    labels = [
        "R1 Thin",
        "R2 Thick",
        "R3 Double-thick",
        "Full TMS",
    ]

    full_vals = [
        *(
            s["Brep_T"]
            for s in full_sections
        ),
        full_all["Brep_T"],
    ]

    central_vals = [
        *(
            s["Brep_T"]
            for s in central_sections
        ),
        central_all["Brep_T"],
    ]

    x = np.arange(
        len(labels),
        dtype=float,
    )

    width = 0.35

    fig, ax = plt.subplots(
        figsize=(11.0, 6.3)
    )

    bars1 = ax.bar(
        x - width / 2,
        full_vals,
        width,
        label="Full steel area",
    )

    bars2 = ax.bar(
        x + width / 2,
        central_vals,
        width,
        label=(
            r"Central steel: "
            r"$10\leq|x|\leq3730$ mm, "
            r"$|y|\leq1780$ mm"
        ),
    )

    for bars in (
        bars1,
        bars2,
    ):
        for bar in bars:
            value = bar.get_height()

            ax.text(
                bar.get_x()
                + bar.get_width() / 2,
                value,
                f"{value:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    ax.set_ylabel(
        "Thickness-weighted representative |B| [T]"
    )

    ax.set_title(
        "Full Steel vs Central Steel Representative Field\n"
        "Same plate-thickness geometry and weighting; only transverse selection changes"
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.22,
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def write_summary(
    path,
    central_sections,
    central_all,
    full_sections,
    full_all,
    data_rows,
    full_rows,
    central_rows,
):
    lines = [
        "TMS CENTRAL-STEEL REPRESENTATIVE-FIELD ANALYSIS",
        "=" * 48,
        "",
        "MAIN ANALYSIS CHAIN",
        "  Full steel plate -> Central steel region -> Flat central region",
        "",
        "THIS SCRIPT",
        "  Main result #2: central steel region",
        "",
        "CENTRAL TRANSVERSE REGION",
        "  10 <= |x| <= 3730 mm",
        "  |y| <= 1780 mm",
        "  20 mm centre gap (-10 < x < +10 mm) excluded",
        "",
        "SELECTION",
        "  Geometry only",
        "  No |B| threshold",
        "  No smoothing",
        "",
        "LONGITUDINAL PLATE-THICKNESS REGIONS",
        "  R1 Thin         : 34 x 15 mm",
        "  R2 Thick        : 22 x 40 mm",
        "  R3 Double-thick : 24 x 80 mm",
        "  Total steel thickness = 3310 mm",
        "",
        "METHOD",
        "  For each physical plate:",
        "      B_i = mean |B| over selected central-steel cells",
        "",
        "  Representative value:",
        "      B_rep = sum_i(B_i * t_i) / sum_i(t_i)",
        "",
        "CENTRAL-STEEL RESULTS",
    ]

    for sec in central_sections:
        lines.append(
            f"  {sec['section']:<18} "
            f"B_rep={sec['Brep_T']:.6f} T  "
            f"n={sec['n_plates']}  "
            f"steel={sec['total_steel_thickness_mm']:.0f} mm"
        )

    lines.append(
        f"  {'Full TMS':<18} "
        f"B_rep={central_all['Brep_T']:.6f} T"
    )

    lines.extend(
        [
            "",
            "FULL-STEEL COMPARISON",
        ]
    )

    for sec in full_sections:
        central_sec = next(
            s
            for s in central_sections
            if s["section"] == sec["section"]
        )

        drop = (
            100.0
            * (
                sec["Brep_T"]
                - central_sec["Brep_T"]
            )
            / sec["Brep_T"]
        )

        lines.append(
            f"  {sec['section']:<18} "
            f"full={sec['Brep_T']:.6f} T  "
            f"central={central_sec['Brep_T']:.6f} T  "
            f"change={drop:.2f}% lower"
        )

    full_drop = (
        100.0
        * (
            full_all["Brep_T"]
            - central_all["Brep_T"]
        )
        / full_all["Brep_T"]
    )

    lines.append(
        f"  {'Full TMS':<18} "
        f"full={full_all['Brep_T']:.6f} T  "
        f"central={central_all['Brep_T']:.6f} T  "
        f"change={full_drop:.2f}% lower"
    )

    lines.extend(
        [
            "",
            "ROW COUNTS",
            f"  Mapper data rows: {data_rows:,}",
            f"  Full-steel plate rows: {full_rows:,}",
            f"  Central-steel plate rows: {central_rows:,}",
            "",
            "IMPORTANT CAVEAT",
            "  The 80 mm source-map offsets use the supplied dataset's 135 mm pitch.",
            "  The current simulation geometry has been observed to use 130 mm,",
            "  so that placement mismatch remains a separate geometry issue.",
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
            "Canonical TMS central-steel representative-field analysis."
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
        default="plots/TMS_central_steel_representative_B",
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
        "TMS CENTRAL-STEEL REPRESENTATIVE-FIELD ANALYSIS"
    )

    print(
        "==============================================="
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
        "Longitudinal plate-thickness regions:"
    )

    print(
        "  R1 / Thin         : 34 x 15 mm"
    )

    print(
        "  R2 / Thick        : 22 x 40 mm"
    )

    print(
        "  R3 / Double-thick : 24 x 80 mm"
    )

    print(
        "  Total steel thickness = 3310 mm"
    )

    print()

    print(
        "Primary definition:"
    )

    print(
        "  same plate geometry and thickness weighting as full-steel analysis"
    )

    print(
        "  only transverse selection changes"
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
        full_acc,
        central_acc,
        data_rows,
        full_rows,
        central_rows,
    ) = stream_mapper(
        mapper_path,
        shape,
        header,
        plates,
    )

    full_plate_rows = plate_rows(
        plates,
        full_acc,
    )

    central_plate_rows = plate_rows(
        plates,
        central_acc,
    )

    full_sections, full_all = summarize(
        full_plate_rows
    )

    central_sections, central_all = summarize(
        central_plate_rows
    )

    write_csv(
        outdir / "central_steel_plate_results.csv",
        central_plate_rows,
    )

    write_csv(
        outdir / "central_steel_section_results.csv",
        central_sections + [central_all],
    )

    write_csv(
        outdir / "full_steel_comparison_plate_results.csv",
        full_plate_rows,
    )

    write_csv(
        outdir / "full_steel_comparison_section_results.csv",
        full_sections + [full_all],
    )

    plot_region_definition(
        outdir / "01_central_steel_region_definition.png"
    )

    plot_central_profile(
        outdir / "02_central_steel_plate_profile.png",
        central_plate_rows,
        central_sections,
        central_all,
    )

    plot_comparison(
        outdir / "03_full_vs_central_steel_comparison.png",
        full_sections,
        central_sections,
        full_all,
        central_all,
    )

    write_summary(
        outdir / "central_steel_representative_B_summary.txt",
        central_sections,
        central_all,
        full_sections,
        full_all,
        data_rows,
        full_rows,
        central_rows,
    )

    print()

    print(
        "CENTRAL-STEEL MAIN RESULTS"
    )

    print(
        "--------------------------"
    )

    for sec in central_sections:
        print(
            f"{sec['section']:<18} "
            f"B_rep={sec['Brep_T']:.4f} T  "
            f"n={sec['n_plates']:2d}  "
            f"steel={sec['total_steel_thickness_mm']:.0f} mm"
        )

    print()

    print(
        "Central-steel thickness-weighted "
        f"B_rep = {central_all['Brep_T']:.4f} T"
    )

    print()

    print(
        "FULL STEEL -> CENTRAL STEEL"
    )

    print(
        "---------------------------"
    )

    for full_sec in full_sections:
        central_sec = next(
            s
            for s in central_sections
            if s["section"] == full_sec["section"]
        )

        change = (
            100.0
            * (
                full_sec["Brep_T"]
                - central_sec["Brep_T"]
            )
            / full_sec["Brep_T"]
        )

        print(
            f"{full_sec['section']:<18} "
            f"{full_sec['Brep_T']:.6f} -> "
            f"{central_sec['Brep_T']:.6f} T "
            f"({change:.2f}% lower)"
        )

    full_change = (
        100.0
        * (
            full_all["Brep_T"]
            - central_all["Brep_T"]
        )
        / full_all["Brep_T"]
    )

    print(
        f"{'Full TMS':<18} "
        f"{full_all['Brep_T']:.6f} -> "
        f"{central_all['Brep_T']:.6f} T "
        f"({full_change:.2f}% lower)"
    )

    print()

    print(
        "Saved:"
    )

    print(
        f"  {outdir / '01_central_steel_region_definition.png'}"
    )

    print(
        f"  {outdir / '02_central_steel_plate_profile.png'}"
    )

    print(
        f"  {outdir / '03_full_vs_central_steel_comparison.png'}"
    )

    print(
        f"  {outdir / 'central_steel_plate_results.csv'}"
    )

    print(
        f"  {outdir / 'central_steel_section_results.csv'}"
    )

    print(
        f"  {outdir / 'central_steel_representative_B_summary.txt'}"
    )


if __name__ == "__main__":
    main()
