#!/usr/bin/env python3
"""
TMS_FULL_STEEL_REPRESENTATIVE_B.py

Canonical MAIN representative-field analysis for the DUNE TMS full steel plate.

PHYSICS DEFINITION
------------------
The TMS is divided by the physical steel-plate thickness geometry:

    R1 / Thin         : 34 plates x 15 mm
    R2 / Thick        : 22 plates x 40 mm
    R3 / Double-thick : 24 plates x 80 mm

For each physical plate:
  1. Select mapper cells using the full-steel transverse geometry mask.
  2. Restrict z to the physical thickness span of that plate.
  3. Compute the per-plate mean |B| with NO magnetic-field threshold.
  4. Compute section representative values using physical plate-thickness
     weighting.
  5. Combine all 80 physical plates into the full-steel representative field:

         B_rep = sum_i( <|B|>_i * t_i ) / sum_i(t_i)

This script defines the FULL-STEEL baseline result in the analysis chain:

    Full steel plate  ->  Central steel region  ->  Flat central region

The separate physical-plate automatic-selection analysis and Maxwell/AEDT
comparison are independent CROSS-CHECKS; they do not replace this definition.

IMPORTANT
---------
- The primary result is geometry-selected.
- No |B| threshold is used in the primary result.
- No smoothing is used in the numerical calculation.
- Threshold variants (>0 and >0.01 T) are diagnostics only.
- The current transverse steel mask is still a first-pass geometry mask until
  exact CAD/GDML shoulder/opening dimensions are fully confirmed.
- The source field uses 135 mm pitch in the 80 mm family, whereas the current
  simulation GDML has been observed to use 130 mm. This remains a geometry
  consistency issue for final production interpretation.

Expected historical reference from the current mapper/mask configuration:
    Thin         ~ 1.447 T
    Thick        ~ 1.366 T
    Double-thick ~ 1.357 T
    Full steel   ~ 1.373 T

The script does NOT hard-code those values; they are recomputed from Mapper.txt.
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

SECTION_ORDER = ("Thin", "Thick", "Double-thick")
SECTION_TO_REGION = {
    "Thin": "R1",
    "Thick": "R2",
    "Double-thick": "R3",
}
SECTION_MARKERS = {
    "Thin": "o",
    "Thick": "s",
    "Double-thick": "^",
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

            vals = s.split()

            if len(vals) != 6:
                raise ValueError(
                    "Expected first non-comment line to contain "
                    "offset_x offset_y offset_z dx dy dz"
                )

            header = tuple(float(v) for v in vals)
            break

    if header is None:
        raise ValueError("Could not find Mapper.txt grid header.")

    return shape, header


def make_plates(local_z0_mm: float):
    plates = []
    number = 1

    for section, thickness, starts in (
        ("Thin", 15.0, THIN_STARTS_MM),
        ("Thick", 40.0, THICK_STARTS_MM),
        ("Double-thick", 80.0, DOUBLE_STARTS_MM),
    ):
        for source_start in starts:
            local_start = local_z0_mm + float(source_start)
            local_end = local_start + thickness

            plates.append(
                Plate(
                    number=number,
                    section=section,
                    thickness_mm=thickness,
                    source_start_mm=float(source_start),
                    local_start_mm=local_start,
                    local_end_mm=local_end,
                    local_center_mm=0.5 * (local_start + local_end),
                )
            )

            number += 1

    return plates


def verify_plate_definition(plates):
    groups = {
        section: [p for p in plates if p.section == section]
        for section in SECTION_ORDER
    }

    assert len(groups["Thin"]) == 34
    assert len(groups["Thick"]) == 22
    assert len(groups["Double-thick"]) == 24
    assert len(plates) == 80

    assert math.isclose(groups["Thin"][0].source_start_mm, 0.0)
    assert math.isclose(groups["Thin"][-1].source_start_mm, 2145.0)

    assert math.isclose(groups["Thick"][0].source_start_mm, 2210.0)
    assert math.isclose(groups["Thick"][-1].source_start_mm, 4100.0)

    assert math.isclose(groups["Double-thick"][0].source_start_mm, 4190.0)
    assert math.isclose(groups["Double-thick"][-1].source_start_mm, 7295.0)

    total_steel_mm = sum(p.thickness_mm for p in plates)
    assert math.isclose(total_steel_mm, 3310.0)


def steel_xy_mask(x: float, y: float, args) -> bool:
    """
    First-pass reflected transverse steel geometry mask.

    The same primary mask used in the previous grand full-steel analysis is
    retained here so the numerical definition remains consistent.
    """
    qx = abs(x)
    qy = abs(y)

    if qx < args.x_min or qx > args.x_max:
        return False

    if qy > args.y_max:
        return False

    # Upper/lower shoulder region.
    if qy > args.shoulder_y:
        if qx < args.shoulder_x_min or qx > args.shoulder_x_max:
            return False

    # Opening/cut-out.
    if (
        args.opening_x_min <= qx <= args.opening_x_max
        and args.opening_y_min <= qy <= args.opening_y_max
    ):
        return False

    return True


def build_z_to_plate(plates, z_values, tolerance=1e-6):
    """
    Assign each regular mapper z-plane to a physical steel plate.

    z planes in air gaps map to None.
    """
    mapping = {}

    for z in z_values:
        plate_index = None

        for i, plate in enumerate(plates):
            if (
                z >= plate.local_start_mm - tolerance
                and z <= plate.local_end_mm + tolerance
            ):
                plate_index = i
                break

        mapping[round(float(z), 6)] = plate_index

    return mapping


# ---------------------------------------------------------------------------
# Streaming accumulation
# ---------------------------------------------------------------------------

def new_accumulator(nplates):
    return {
        "count_geom": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bperp_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bx_geom": np.zeros(nplates, dtype=np.float64),
        "sum_by_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bz_geom": np.zeros(nplates, dtype=np.float64),

        # Diagnostics only:
        "count_gt0": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_gt0": np.zeros(nplates, dtype=np.float64),
        "count_gt001": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_gt001": np.zeros(nplates, dtype=np.float64),
        "zero_count": np.zeros(nplates, dtype=np.int64),
    }


def accumulate_one(acc, idx, bx, by, bz, bmag, zero_tolerance):
    bperp = math.sqrt(bx * bx + by * by)

    acc["count_geom"][idx] += 1
    acc["sum_bmag_geom"][idx] += bmag
    acc["sum_bperp_geom"][idx] += bperp
    acc["sum_bx_geom"][idx] += bx
    acc["sum_by_geom"][idx] += by
    acc["sum_bz_geom"][idx] += bz

    if bmag <= zero_tolerance:
        acc["zero_count"][idx] += 1

    if bmag > zero_tolerance:
        acc["count_gt0"][idx] += 1
        acc["sum_bmag_gt0"][idx] += bmag

    if bmag > 0.01:
        acc["count_gt001"][idx] += 1
        acc["sum_bmag_gt001"][idx] += bmag


def stream_mapper(mapper_path: Path, plates, shape, header, args):
    _, _, oz, _, _, dz = header

    if shape is None:
        max_end = max(p.local_end_mm for p in plates)
        nz = int(math.ceil((max_end - oz) / dz)) + 1
    else:
        _, _, nz = shape

    z_values = oz + np.arange(nz, dtype=float) * dz
    z_to_plate = build_z_to_plate(plates, z_values)

    acc = new_accumulator(len(plates))

    header_seen = False
    data_rows = 0
    xy_steel_rows = 0
    plate_rows = 0
    gap_rows = 0

    with mapper_path.open("r") as f:
        for line in f:
            s = line.strip()

            if not s or s.startswith("#"):
                continue

            vals = s.split()

            if not header_seen:
                if len(vals) != 6:
                    raise ValueError("Expected six-number mapper grid header.")
                header_seen = True
                continue

            if len(vals) < 6:
                continue

            x, y, z, bx, by, bz = map(float, vals[:6])

            if len(vals) >= 7:
                bmag = float(vals[6])
            else:
                bmag = math.sqrt(bx * bx + by * by + bz * bz)

            data_rows += 1

            if not steel_xy_mask(x, y, args):
                continue

            xy_steel_rows += 1

            plate_idx = z_to_plate.get(round(z, 6))

            if plate_idx is None:
                gap_rows += 1
                continue

            plate_rows += 1

            accumulate_one(
                acc,
                plate_idx,
                bx,
                by,
                bz,
                bmag,
                args.zero_tolerance,
            )

    return {
        "acc": acc,
        "data_rows": data_rows,
        "steel_rows": xy_steel_rows,
        "plate_rows": plate_rows,
        "gap_rows": gap_rows,
    }


# ---------------------------------------------------------------------------
# Plate / section summaries
# ---------------------------------------------------------------------------

def safe_mean(sum_arr, count_arr, i):
    n = int(count_arr[i])

    if n == 0:
        return float("nan")

    return float(sum_arr[i] / n)


def build_plate_rows(plates, acc):
    rows = []

    for i, plate in enumerate(plates):
        n_geom = int(acc["count_geom"][i])

        mean_bmag = safe_mean(
            acc["sum_bmag_geom"],
            acc["count_geom"],
            i,
        )

        mean_bperp = safe_mean(
            acc["sum_bperp_geom"],
            acc["count_geom"],
            i,
        )

        mean_bx = safe_mean(
            acc["sum_bx_geom"],
            acc["count_geom"],
            i,
        )

        mean_by = safe_mean(
            acc["sum_by_geom"],
            acc["count_geom"],
            i,
        )

        mean_bz = safe_mean(
            acc["sum_bz_geom"],
            acc["count_geom"],
            i,
        )

        mean_gt0 = safe_mean(
            acc["sum_bmag_gt0"],
            acc["count_gt0"],
            i,
        )

        mean_gt001 = safe_mean(
            acc["sum_bmag_gt001"],
            acc["count_gt001"],
            i,
        )

        zero_fraction = (
            float(acc["zero_count"][i] / n_geom)
            if n_geom > 0
            else float("nan")
        )

        thickness_m = plate.thickness_mm / 1000.0

        rows.append(
            {
                "plate": plate.number,
                "region": SECTION_TO_REGION[plate.section],
                "section": plate.section,
                "thickness_mm": plate.thickness_mm,
                "source_start_mm": plate.source_start_mm,
                "local_start_mm": plate.local_start_mm,
                "local_end_mm": plate.local_end_mm,
                "local_center_mm": plate.local_center_mm,
                "n_geometry_cells": n_geom,
                "zero_fraction": zero_fraction,
                "mean_Bmag_geometry_T": mean_bmag,
                "mean_Bperp_geometry_T": mean_bperp,
                "mean_Bx_geometry_T": mean_bx,
                "mean_By_geometry_T": mean_by,
                "mean_Bz_geometry_T": mean_bz,
                "mean_Bmag_gt0_T": mean_gt0,
                "mean_Bmag_gt001_T": mean_gt001,
                "Bmag_times_thickness_Tm": mean_bmag * thickness_m,
                "Bperp_times_thickness_Tm": mean_bperp * thickness_m,
            }
        )

    return rows


def weighted_mean(rows, field):
    valid = [
        r for r in rows
        if np.isfinite(r[field])
    ]

    denominator = sum(
        r["thickness_mm"]
        for r in valid
    )

    if denominator == 0:
        return float("nan")

    numerator = sum(
        r[field] * r["thickness_mm"]
        for r in valid
    )

    return numerator / denominator


def summarize_section(rows, section):
    selected = [
        r for r in rows
        if r["section"] == section
        and np.isfinite(r["mean_Bmag_geometry_T"])
    ]

    if not selected:
        raise ValueError(f"No valid plates for section {section}.")

    total_mm = sum(
        r["thickness_mm"]
        for r in selected
    )

    return {
        "region": SECTION_TO_REGION[section],
        "section": section,
        "thickness_mm": selected[0]["thickness_mm"],
        "n_plates": len(selected),
        "source_z_start_mm": min(
            r["source_start_mm"]
            for r in selected
        ),
        "source_z_last_start_mm": max(
            r["source_start_mm"]
            for r in selected
        ),
        "local_z_start_mm": min(
            r["local_start_mm"]
            for r in selected
        ),
        "local_z_end_mm": max(
            r["local_end_mm"]
            for r in selected
        ),
        "total_steel_thickness_mm": total_mm,
        "total_steel_thickness_m": total_mm / 1000.0,

        # PRIMARY
        "Brep_geometry_T": weighted_mean(
            selected,
            "mean_Bmag_geometry_T",
        ),
        "Bperp_rep_geometry_T": weighted_mean(
            selected,
            "mean_Bperp_geometry_T",
        ),

        # DIAGNOSTIC ONLY
        "Brep_gt0_T": weighted_mean(
            selected,
            "mean_Bmag_gt0_T",
        ),
        "Brep_gt001_T": weighted_mean(
            selected,
            "mean_Bmag_gt001_T",
        ),

        "integral_absB_Tm": sum(
            r["Bmag_times_thickness_Tm"]
            for r in selected
        ),
        "integral_Bperp_Tm": sum(
            r["Bperp_times_thickness_Tm"]
            for r in selected
        ),
        "mean_zero_fraction": float(
            np.nanmean([
                r["zero_fraction"]
                for r in selected
            ])
        ),
    }


def summarize_sections(rows):
    return [
        summarize_section(rows, section)
        for section in SECTION_ORDER
    ]


def summarize_full(rows):
    selected = [
        r for r in rows
        if np.isfinite(r["mean_Bmag_geometry_T"])
    ]

    total_mm = sum(
        r["thickness_mm"]
        for r in selected
    )

    return {
        "region": "Full",
        "section": "Full steel",
        "n_plates": len(selected),
        "total_steel_thickness_mm": total_mm,
        "total_steel_thickness_m": total_mm / 1000.0,

        # PRIMARY full-steel representative field
        "Brep_geometry_T": weighted_mean(
            selected,
            "mean_Bmag_geometry_T",
        ),
        "Bperp_rep_geometry_T": weighted_mean(
            selected,
            "mean_Bperp_geometry_T",
        ),

        # DIAGNOSTIC ONLY
        "Brep_gt0_T": weighted_mean(
            selected,
            "mean_Bmag_gt0_T",
        ),
        "Brep_gt001_T": weighted_mean(
            selected,
            "mean_Bmag_gt001_T",
        ),

        "integral_absB_Tm": sum(
            r["Bmag_times_thickness_Tm"]
            for r in selected
        ),
        "integral_Bperp_Tm": sum(
            r["Bperp_times_thickness_Tm"]
            for r in selected
        ),
    }


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(path, rows):
    if not rows:
        return

    fieldnames = list(rows[0].keys())

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_full_steel_profile(path, plate_rows, sections, full):
    fig, ax = plt.subplots(figsize=(13.5, 7.0))

    for sec in sections:
        section = sec["section"]

        selected = [
            r for r in plate_rows
            if r["section"] == section
        ]

        z = np.asarray(
            [
                r["local_center_mm"]
                for r in selected
            ],
            dtype=float,
        )

        b = np.asarray(
            [
                r["mean_Bmag_geometry_T"]
                for r in selected
            ],
            dtype=float,
        )

        ax.plot(
            z,
            b,
            marker=SECTION_MARKERS[section],
            markersize=5,
            linewidth=1.4,
            label=(
                f"{sec['region']} — {section} "
                f"({sec['thickness_mm']:.0f} mm): "
                f"B_rep={sec['Brep_geometry_T']:.3f} T"
            ),
        )

        ax.hlines(
            sec["Brep_geometry_T"],
            sec["local_z_start_mm"],
            sec["local_z_end_mm"],
            linestyles="--",
            linewidth=1.4,
        )

    thick_start = next(
        r["local_start_mm"]
        for r in plate_rows
        if r["section"] == "Thick"
    )

    double_start = next(
        r["local_start_mm"]
        for r in plate_rows
        if r["section"] == "Double-thick"
    )

    ax.axvline(
        thick_start,
        linestyle=":",
        linewidth=1.2,
    )

    ax.axvline(
        double_start,
        linestyle=":",
        linewidth=1.2,
    )

    ax.axhline(
        full["Brep_geometry_T"],
        linestyle="-.",
        linewidth=2.0,
        label=(
            "Full-steel thickness-weighted "
            f"B_rep = {full['Brep_geometry_T']:.3f} T"
        ),
    )

    ax.set_title(
        "TMS Representative Magnetic Field — Full Steel Plate\n"
        "Regions defined by 15 mm, 40 mm, and 80 mm plate-thickness geometry",
        fontsize=15,
    )

    ax.set_xlabel(
        "TMS local Z [mm]",
        fontsize=13,
    )

    ax.set_ylabel(
        "Per-plate mean |B| [T]",
        fontsize=13,
    )

    ax.text(
        0.01,
        0.02,
        (
            "Primary full-steel result: geometry selection, "
            "no B threshold, no smoothing"
        ),
        transform=ax.transAxes,
        fontsize=9,
    )

    ax.tick_params(
        axis="both",
        labelsize=11,
    )

    ax.grid(
        True,
        alpha=0.25,
    )

    ax.legend(
        loc="best",
        fontsize=9,
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_bending_power(path, sections):
    labels = [
        (
            f"{s['region']}\n"
            f"{s['section']}\n"
            f"({s['thickness_mm']:.0f} mm)"
        )
        for s in sections
    ]

    values = [
        s["integral_Bperp_Tm"]
        for s in sections
    ]

    fig, ax = plt.subplots(
        figsize=(9.6, 6.0)
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
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.3f} T·m",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax.set_title(
        "Full-Steel Bending-Power Proxy by Plate-Thickness Region",
        fontsize=15,
    )

    ax.set_ylabel(
        r"$\sum \langle B_{\perp}\rangle_i\,t_i$ [T·m]",
        fontsize=12,
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


def plot_threshold_diagnostic(path, sections, full):
    rows = sections + [full]

    labels = [
        (
            r["region"]
            if r["region"] != "Full"
            else "Full steel"
        )
        for r in rows
    ]

    geometry = [
        r["Brep_geometry_T"]
        for r in rows
    ]

    gt0 = [
        r["Brep_gt0_T"]
        for r in rows
    ]

    gt001 = [
        r["Brep_gt001_T"]
        for r in rows
    ]

    x = np.arange(
        len(labels),
        dtype=float,
    )

    width = 0.25

    fig, ax = plt.subplots(
        figsize=(11.0, 6.2)
    )

    ax.bar(
        x - width,
        geometry,
        width,
        label="Geometry selected; no threshold",
    )

    ax.bar(
        x,
        gt0,
        width,
        label="|B| > 0",
    )

    ax.bar(
        x + width,
        gt001,
        width,
        label="|B| > 0.01 T",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    ax.set_ylabel(
        "Representative |B| [T]"
    )

    ax.set_title(
        "Threshold Sensitivity Diagnostic\n"
        "Primary result is geometry-selected; threshold variants are cross-checks only"
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    ax.legend(
        fontsize=9
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def write_summary(path, args, sections, full, stream_stats):
    lines = []

    lines.append(
        "TMS FULL-STEEL REPRESENTATIVE-FIELD ANALYSIS"
    )

    lines.append(
        "=" * 48
    )

    lines.append("")

    lines.append(
        "PRIMARY RESULT IN THE REPRESENTATIVE-FIELD ANALYSIS CHAIN"
    )

    lines.append(
        "  Full steel plate -> Central steel region -> Flat central region"
    )

    lines.append("")

    lines.append(
        "PHYSICAL PLATE-THICKNESS REGIONS"
    )

    lines.append(
        "  R1 / Thin         : 34 x 15 mm, starts 0 ... 2145 mm"
    )

    lines.append(
        "  R2 / Thick        : 22 x 40 mm, starts 2210 ... 4100 mm"
    )

    lines.append(
        "  R3 / Double-thick : 24 x 80 mm, starts 4190 ... 7295 mm"
    )

    lines.append(
        "  Total physical steel thickness = 3310 mm"
    )

    lines.append("")

    lines.append(
        "PRIMARY DEFINITION"
    )

    lines.append(
        "  Geometry determines steel cells; no B threshold."
    )

    lines.append(
        "  One mean |B| is calculated for each physical plate."
    )

    lines.append(
        "  Plate/section/full values are weighted by physical steel thickness."
    )

    lines.append(
        "  No smoothing is used."
    )

    lines.append("")

    lines.append(
        "FORMULA"
    )

    lines.append(
        "  B_rep = sum_i( <|B|>_i * t_i ) / sum_i(t_i)"
    )

    lines.append("")

    lines.append(
        "RESULTS"
    )

    for sec in sections:
        lines.append(
            f"  {sec['region']} / {sec['section']:<13} "
            f"B_rep={sec['Brep_geometry_T']:.6f} T  "
            f"n={sec['n_plates']}  "
            f"steel={sec['total_steel_thickness_mm']:.0f} mm"
        )

    lines.append("")

    lines.append(
        f"  Full steel: B_rep={full['Brep_geometry_T']:.6f} T"
    )

    lines.append(
        f"  Full steel thickness={full['total_steel_thickness_mm']:.0f} mm"
    )

    lines.append("")

    lines.append(
        "BENDING-POWER PROXY"
    )

    for sec in sections:
        lines.append(
            f"  {sec['region']} / {sec['section']:<13} "
            f"Integral(B_perp dl)~{sec['integral_Bperp_Tm']:.6f} T*m"
        )

    lines.append(
        f"  Full steel Integral(B_perp dl)~{full['integral_Bperp_Tm']:.6f} T*m"
    )

    lines.append("")

    lines.append(
        "THRESHOLD DIAGNOSTICS — NOT THE PRIMARY DEFINITION"
    )

    for sec in sections:
        lines.append(
            f"  {sec['region']} / {sec['section']:<13} "
            f"geometry={sec['Brep_geometry_T']:.6f} T, "
            f">0={sec['Brep_gt0_T']:.6f} T, "
            f">0.01={sec['Brep_gt001_T']:.6f} T, "
            f"mean zero fraction={100.0 * sec['mean_zero_fraction']:.2f}%"
        )

    lines.append("")

    lines.append(
        "STREAMING DIAGNOSTICS"
    )

    lines.append(
        f"  Mapper data rows: {stream_stats['data_rows']:,}"
    )

    lines.append(
        f"  XY steel-mask rows: {stream_stats['steel_rows']:,}"
    )

    lines.append(
        f"  Rows inside physical plate z spans: {stream_stats['plate_rows']:,}"
    )

    lines.append(
        f"  Rows in z gaps after XY mask: {stream_stats['gap_rows']:,}"
    )

    lines.append("")

    lines.append(
        "IMPORTANT CAVEATS"
    )

    lines.append(
        "  * The transverse steel mask is still first-pass until exact "
        "GDML/CAD shoulder/opening dimensions are confirmed."
    )

    lines.append(
        "  * The current field dataset uses 135 mm pitch in the 80 mm section; "
        "the current simulation GDML has been observed to use 130 mm."
    )

    lines.append(
        "  * Maxwell/AEDT line scans and the automatic physical-plate analysis "
        "are independent cross-checks, not replacements for this definition."
    )

    lines.append(
        "  * A true muon bend requires track-dependent integral |B x u_track| dl."
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description=(
            "Canonical TMS full-steel representative-field analysis: "
            "physical plate-thickness regions, geometry selection, "
            "no B threshold, no smoothing."
        )
    )

    p.add_argument(
        "mapper",
        nargs="?",
        default="Mapper.txt",
    )

    p.add_argument(
        "-o",
        "--output-dir",
        default="plots/TMS_full_steel_representative_B",
    )

    p.add_argument(
        "--zero-tolerance",
        type=float,
        default=1e-12,
        help="Diagnostic definition of exact/near-zero |B|.",
    )

    # Same first-pass transverse steel mask as the previous grand analysis.
    p.add_argument(
        "--x-min",
        type=float,
        default=10.0,
    )

    p.add_argument(
        "--x-max",
        type=float,
        default=3730.0,
    )

    p.add_argument(
        "--y-max",
        type=float,
        default=2350.0,
    )

    p.add_argument(
        "--shoulder-y",
        type=float,
        default=1780.0,
    )

    p.add_argument(
        "--shoulder-x-min",
        type=float,
        default=250.0,
    )

    p.add_argument(
        "--shoulder-x-max",
        type=float,
        default=3480.0,
    )

    p.add_argument(
        "--opening-x-min",
        type=float,
        default=1620.0,
    )

    p.add_argument(
        "--opening-x-max",
        type=float,
        default=2120.0,
    )

    p.add_argument(
        "--opening-y-min",
        type=float,
        default=1780.0,
    )

    p.add_argument(
        "--opening-y-max",
        type=float,
        default=2000.0,
    )

    return p


def main():
    args = build_parser().parse_args()

    mapper_path = Path(args.mapper)
    outdir = Path(args.output_dir)

    if not mapper_path.exists():
        raise FileNotFoundError(mapper_path)

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    shape, header = read_mapper_header(
        mapper_path
    )

    _, _, oz, _, _, _ = header

    plates = make_plates(
        local_z0_mm=oz
    )

    verify_plate_definition(
        plates
    )

    print(
        "TMS FULL-STEEL REPRESENTATIVE-FIELD ANALYSIS"
    )

    print(
        "============================================"
    )

    print()

    print(
        "Physical plate-thickness regions:"
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
        "  geometry-selected full steel area"
    )

    print(
        "  no B threshold"
    )

    print(
        "  no smoothing"
    )

    print(
        "  per-plate mean |B| -> physical thickness weighting"
    )

    print()

    print(
        f"Mapper local source origin z0 = {oz:g} mm"
    )

    print()

    print(
        "Reading Mapper.txt in one pass..."
    )

    stats = stream_mapper(
        mapper_path,
        plates,
        shape,
        header,
        args,
    )

    plate_rows = build_plate_rows(
        plates,
        stats["acc"],
    )

    sections = summarize_sections(
        plate_rows
    )

    full = summarize_full(
        plate_rows
    )

    write_csv(
        outdir / "full_plate_results.csv",
        plate_rows,
    )

    write_csv(
        outdir / "full_section_results.csv",
        sections + [full],
    )

    plot_full_steel_profile(
        outdir / "01_full_steel_plate_profile.png",
        plate_rows,
        sections,
        full,
    )

    plot_bending_power(
        outdir / "02_full_steel_bending_power_by_section.png",
        sections,
    )

    plot_threshold_diagnostic(
        outdir / "03_threshold_sensitivity_diagnostic.png",
        sections,
        full,
    )

    write_summary(
        outdir / "full_steel_representative_B_summary.txt",
        args,
        sections,
        full,
        stats,
    )

    print()

    print(
        "FULL-STEEL MAIN RESULTS"
    )

    print(
        "-----------------------"
    )

    for sec in sections:
        print(
            f"{sec['region']} / {sec['section']:<13} "
            f"B_rep={sec['Brep_geometry_T']:.4f} T  "
            f"n={sec['n_plates']:2d}  "
            f"steel={sec['total_steel_thickness_mm']:.0f} mm"
        )

    print()

    print(
        "Full-steel thickness-weighted "
        f"B_rep = {full['Brep_geometry_T']:.4f} T"
    )

    print(
        "Full-steel integral(B_perp dl) proxy = "
        f"{full['integral_Bperp_Tm']:.4f} T*m"
    )

    print()

    print(
        "Threshold variants are diagnostic only."
    )

    print()

    print(
        "Saved:"
    )

    print(
        f"  {outdir / '01_full_steel_plate_profile.png'}"
    )

    print(
        f"  {outdir / '02_full_steel_bending_power_by_section.png'}"
    )

    print(
        f"  {outdir / '03_threshold_sensitivity_diagnostic.png'}"
    )

    print(
        f"  {outdir / 'full_plate_results.csv'}"
    )

    print(
        f"  {outdir / 'full_section_results.csv'}"
    )

    print(
        f"  {outdir / 'full_steel_representative_B_summary.txt'}"
    )


if __name__ == "__main__":
    main()
