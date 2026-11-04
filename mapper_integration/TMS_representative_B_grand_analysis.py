#!/usr/bin/env python3
"""
Grand TMS representative-field analysis
=======================================

This script implements the corrected analysis logic discussed after the expert
feedback.

PRIMARY PHYSICS CHANGES
-----------------------
1. Divide the TMS by the REAL plate-thickness sections:
     Thin         : 34 plates x 15 mm, source starts 0 ... 2145 mm
     Thick        : 22 plates x 40 mm, source starts 2210 ... 4100 mm
     Double-thick : 24 plates x 80 mm, source starts 4190 ... 7295 mm

2. Calculate one representative field for EACH physical plate.

3. Calculate each thickness section independently:
       B_section = sum_i(B_i * t_i) / sum_i(t_i)

   Since all plates inside a given section have the same t_i, this reduces
   mathematically to the mean of the plate-level B_i values in that section.

4. Only AFTER the three section values are understood, optionally combine all
   plates into one full-TMS steel-thickness-weighted value:
       B_TMS = sum_i(B_i * t_i) / sum_i(t_i)

5. PRIMARY selection is GEOMETRY-BASED:
       transverse steel mask + physical z thickness of each steel plate.

   No B threshold is used in the primary result.

6. For transparency, ALSO calculate diagnostic versions using:
       |B| > 0
       |B| > 0.01 T

   These are diagnostics only, not the primary definition of steel.

7. NO smoothing is used anywhere in the numerical calculation.

8. Calculate a bending-power proxy for each thickness section:
       I_perp = sum_i(<B_perp>_i * t_i)
   where for a beam-parallel (+z) reference direction:
       B_perp = sqrt(Bx^2 + By^2)

   I_perp has units T*m and is the quantity most directly related to magnetic
   bending strength.  For an actual muon track, the rigorous quantity is
       integral |B x u_track| dl
   along that track.

9. Optional CENTRAL/FIDUCIAL analysis:
   If central x/y bounds are supplied, all calculations are repeated using only
   that central region.  This lets you compare:
       full-steel representative field
       vs
       central/fiducial representative field

   IMPORTANT: Do not choose the central box to force the result toward 1 T.
   Use bounds agreed with the magnet experts / Marco.

OUTPUTS
-------
plots/TMS_representative_B_grand/
  01_full_plate_profile.png
  02_full_bending_power_by_section.png
  03_threshold_sensitivity.png
  04_central_vs_full_sections.png       (only if central region requested)
  05_central_plate_profile.png          (only if central region requested)
  06_central_bending_power_by_section.png (only if central region requested)

  full_plate_results.csv
  full_section_results.csv
  central_plate_results.csv             (if requested)
  central_section_results.csv           (if requested)
  grand_summary.txt

DEPENDENCIES
------------
numpy
matplotlib

No pandas or scipy required.
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

# Source z-coordinate definitions from the current field dataset.
THIN_STARTS_MM = np.arange(34, dtype=float) * 65.0
THICK_STARTS_MM = 2210.0 + np.arange(22, dtype=float) * 90.0
DOUBLE_STARTS_MM = 4190.0 + np.arange(24, dtype=float) * 135.0

SECTION_ORDER = ("Thin", "Thick", "Double-thick")

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
                    "Expected the first non-comment line to be: "
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

    total = sum(p.thickness_mm for p in plates)
    assert math.isclose(total, 3310.0)


def steel_xy_mask(x: float, y: float, args) -> bool:
    """
    Current first-pass reflected transverse steel mask.

    Dimensions explicitly visible in the plate drawing are used where known.
    The shoulder x limits and opening y_max remain first-pass inferred values
    until exact GDML/CAD values are supplied.
    """
    qx = abs(x)
    qy = abs(y)

    if qx < args.x_min or qx > args.x_max:
        return False

    if qy > args.y_max:
        return False

    if qy > args.shoulder_y:
        if qx < args.shoulder_x_min or qx > args.shoulder_x_max:
            return False

    if (
        args.opening_x_min <= qx <= args.opening_x_max
        and args.opening_y_min <= qy <= args.opening_y_max
    ):
        return False

    return True


def central_region_requested(args) -> bool:
    values = (
        args.central_x_min,
        args.central_x_max,
        args.central_y_min,
        args.central_y_max,
    )
    return any(v is not None for v in values)


def validate_central_region(args):
    if not central_region_requested(args):
        return False

    required = {
        "--central-x-min": args.central_x_min,
        "--central-x-max": args.central_x_max,
        "--central-y-max": args.central_y_max,
    }

    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError(
            "Central analysis requested, but missing: "
            + ", ".join(missing)
        )

    if args.central_y_min is None:
        args.central_y_min = 0.0

    if args.central_x_min < 0:
        raise ValueError("--central-x-min must be >= 0.")

    if args.central_x_max <= args.central_x_min:
        raise ValueError(
            "--central-x-max must be greater than --central-x-min."
        )

    if args.central_y_min < 0:
        raise ValueError("--central-y-min must be >= 0.")

    if args.central_y_max <= args.central_y_min:
        raise ValueError(
            "--central-y-max must be greater than --central-y-min."
        )

    return True


def central_xy_mask(x: float, y: float, args) -> bool:
    """
    Symmetric central/fiducial region in |x|, |y|.

    This region is always intersected with the steel geometry mask.
    """
    if not steel_xy_mask(x, y, args):
        return False

    qx = abs(x)
    qy = abs(y)

    return (
        args.central_x_min <= qx <= args.central_x_max
        and args.central_y_min <= qy <= args.central_y_max
    )


def build_z_to_plate(plates, z_values, tolerance=1e-6):
    """
    Assign mapper z planes to physical steel plates.
    Planes in air gaps are mapped to None.
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


def new_accumulator(nplates):
    return {
        "count_geom": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bperp_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bx_geom": np.zeros(nplates, dtype=np.float64),
        "sum_by_geom": np.zeros(nplates, dtype=np.float64),
        "sum_bz_geom": np.zeros(nplates, dtype=np.float64),

        "count_gt0": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_gt0": np.zeros(nplates, dtype=np.float64),

        "count_gt001": np.zeros(nplates, dtype=np.int64),
        "sum_bmag_gt001": np.zeros(nplates, dtype=np.float64),

        "zero_count": np.zeros(nplates, dtype=np.int64),
    }


def accumulate_one(acc, idx, bx, by, bz, bmag, args):
    # Beam-parallel reference direction u = +z.
    # B_perp = |B - (B.u)u| = sqrt(Bx^2 + By^2)
    bperp = math.sqrt(bx * bx + by * by)

    acc["count_geom"][idx] += 1
    acc["sum_bmag_geom"][idx] += bmag
    acc["sum_bperp_geom"][idx] += bperp
    acc["sum_bx_geom"][idx] += bx
    acc["sum_by_geom"][idx] += by
    acc["sum_bz_geom"][idx] += bz

    if bmag <= args.zero_tolerance:
        acc["zero_count"][idx] += 1

    if bmag > args.zero_tolerance:
        acc["count_gt0"][idx] += 1
        acc["sum_bmag_gt0"][idx] += bmag

    if bmag > 0.01:
        acc["count_gt001"][idx] += 1
        acc["sum_bmag_gt001"][idx] += bmag


def stream_mapper(mapper_path: Path, plates, shape, header, args, do_central):
    ox, oy, oz, dx, dy, dz = header

    if shape is None:
        max_end = max(p.local_end_mm for p in plates)
        nz = int(math.ceil((max_end - oz) / dz)) + 1
    else:
        _, _, nz = shape

    z_values = oz + np.arange(nz, dtype=float) * dz
    z_to_plate = build_z_to_plate(plates, z_values)

    full = new_accumulator(len(plates))
    central = new_accumulator(len(plates)) if do_central else None

    header_seen = False
    data_rows = 0
    steel_rows = 0
    plate_rows = 0
    gap_rows = 0
    central_rows = 0

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

            steel_rows += 1

            plate_idx = z_to_plate.get(round(z, 6))
            if plate_idx is None:
                gap_rows += 1
                continue

            plate_rows += 1
            accumulate_one(full, plate_idx, bx, by, bz, bmag, args)

            if do_central and central_xy_mask(x, y, args):
                central_rows += 1
                accumulate_one(
                    central,
                    plate_idx,
                    bx,
                    by,
                    bz,
                    bmag,
                    args,
                )

    return {
        "full": full,
        "central": central,
        "data_rows": data_rows,
        "steel_rows": steel_rows,
        "plate_rows": plate_rows,
        "gap_rows": gap_rows,
        "central_rows": central_rows,
    }


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

        t_m = plate.thickness_mm / 1000.0

        rows.append(
            {
                "plate": plate.number,
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

                "Bmag_times_thickness_Tm": mean_bmag * t_m,
                "Bperp_times_thickness_Tm": mean_bperp * t_m,
            }
        )

    return rows


def weighted_section(rows, section):
    selected = [
        r for r in rows
        if r["section"] == section
        and np.isfinite(r["mean_Bmag_geometry_T"])
    ]

    if not selected:
        raise ValueError(f"No valid plates for section {section}.")

    total_mm = sum(r["thickness_mm"] for r in selected)
    total_m = total_mm / 1000.0

    def weighted(field):
        numerator = sum(
            r[field] * r["thickness_mm"]
            for r in selected
            if np.isfinite(r[field])
        )
        denominator = sum(
            r["thickness_mm"]
            for r in selected
            if np.isfinite(r[field])
        )
        return numerator / denominator if denominator else float("nan")

    bmag_geom = weighted("mean_Bmag_geometry_T")
    bperp_geom = weighted("mean_Bperp_geometry_T")
    bmag_gt0 = weighted("mean_Bmag_gt0_T")
    bmag_gt001 = weighted("mean_Bmag_gt001_T")

    # Scalar bending-power proxy for beam-parallel tracks.
    integral_bperp = sum(
        r["Bperp_times_thickness_Tm"]
        for r in selected
    )

    integral_bmag = sum(
        r["Bmag_times_thickness_Tm"]
        for r in selected
    )

    mean_zero_fraction = float(
        np.nanmean([r["zero_fraction"] for r in selected])
    )

    return {
        "section": section,
        "thickness_mm": selected[0]["thickness_mm"],
        "n_plates": len(selected),
        "source_z_start_mm": min(r["source_start_mm"] for r in selected),
        "source_z_last_start_mm": max(r["source_start_mm"] for r in selected),
        "local_z_start_mm": min(r["local_start_mm"] for r in selected),
        "local_z_end_mm": max(r["local_end_mm"] for r in selected),

        "total_steel_thickness_mm": total_mm,
        "total_steel_thickness_m": total_m,

        "Brep_geometry_T": bmag_geom,
        "Bperp_rep_geometry_T": bperp_geom,

        "Brep_gt0_T": bmag_gt0,
        "Brep_gt001_T": bmag_gt001,

        "integral_absB_Tm": integral_bmag,
        "integral_Bperp_Tm": integral_bperp,

        "mean_zero_fraction": mean_zero_fraction,
    }


def summarize_sections(rows):
    return [
        weighted_section(rows, section)
        for section in SECTION_ORDER
    ]


def summarize_full(rows):
    selected = [
        r for r in rows
        if np.isfinite(r["mean_Bmag_geometry_T"])
    ]

    total_mm = sum(r["thickness_mm"] for r in selected)

    def weighted(field):
        numerator = sum(
            r[field] * r["thickness_mm"]
            for r in selected
            if np.isfinite(r[field])
        )
        denominator = sum(
            r["thickness_mm"]
            for r in selected
            if np.isfinite(r[field])
        )
        return numerator / denominator if denominator else float("nan")

    return {
        "section": "Full TMS",
        "n_plates": len(selected),
        "total_steel_thickness_mm": total_mm,
        "total_steel_thickness_m": total_mm / 1000.0,

        "Brep_geometry_T": weighted("mean_Bmag_geometry_T"),
        "Bperp_rep_geometry_T": weighted("mean_Bperp_geometry_T"),
        "Brep_gt0_T": weighted("mean_Bmag_gt0_T"),
        "Brep_gt001_T": weighted("mean_Bmag_gt001_T"),

        "integral_absB_Tm": sum(
            r["Bmag_times_thickness_Tm"]
            for r in selected
        ),
        "integral_Bperp_Tm": sum(
            r["Bperp_times_thickness_Tm"]
            for r in selected
        ),
    }


def add_bend_angle(summary, momentum_gev):
    if momentum_gev is None:
        summary["approx_bend_angle_rad"] = float("nan")
        summary["approx_bend_angle_deg"] = float("nan")
        return

    # Standard high-energy approximation:
    # theta [rad] ~= 0.299792458 * integral(B_perp dl) [T*m] / p [GeV/c]
    theta = (
        0.299792458
        * summary["integral_Bperp_Tm"]
        / momentum_gev
    )

    summary["approx_bend_angle_rad"] = theta
    summary["approx_bend_angle_deg"] = math.degrees(theta)


def write_csv(path, rows):
    if not rows:
        return

    fieldnames = list(rows[0].keys())

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_plate_profile(path, plate_rows, sections, full, title_prefix):
    fig, ax = plt.subplots(figsize=(13.5, 7.0))

    for sec in sections:
        section = sec["section"]

        selected = [
            r for r in plate_rows
            if r["section"] == section
        ]

        z = np.asarray(
            [r["local_center_mm"] for r in selected],
            dtype=float,
        )

        b = np.asarray(
            [r["mean_Bmag_geometry_T"] for r in selected],
            dtype=float,
        )

        ax.plot(
            z,
            b,
            marker=SECTION_MARKERS[section],
            markersize=5,
            linewidth=1.4,
            label=(
                f"{section} ({sec['thickness_mm']:.0f} mm): "
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

    # Physical section boundaries.
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

    ax.axvline(thick_start, linestyle=":", linewidth=1.2)
    ax.axvline(double_start, linestyle=":", linewidth=1.2)

    ax.axhline(
        full["Brep_geometry_T"],
        linestyle="-.",
        linewidth=2.0,
        label=(
            f"Full TMS thickness-weighted B_rep = "
            f"{full['Brep_geometry_T']:.3f} T"
        ),
    )

    ax.set_title(
        f"{title_prefix}\n"
        "Per-plate mean |B|, then steel-thickness weighting; "
        "geometry selection, no B threshold",
        fontsize=15,
    )

    ax.set_xlabel("TMS local Z [mm]", fontsize=13)
    ax.set_ylabel("Per-plate mean |B| [T]", fontsize=13)
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=9)

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_bending_power(path, sections, title_prefix):
    labels = [
        f"{s['section']}\n({s['thickness_mm']:.0f} mm)"
        for s in sections
    ]

    values = [
        s["integral_Bperp_Tm"]
        for s in sections
    ]

    fig, ax = plt.subplots(figsize=(9.6, 6.0))

    bars = ax.bar(labels, values)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.3f} T·m",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax.set_title(
        f"{title_prefix}\n"
        "Bending-power proxy by thickness section",
        fontsize=15,
    )
    ax.set_ylabel(r"$\sum \langle B_{\perp}\rangle_i\,t_i$  [T·m]", fontsize=12)
    ax.grid(True, axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_threshold_sensitivity(path, sections, full):
    rows = sections + [full]

    labels = [r["section"] for r in rows]
    geometry = [r["Brep_geometry_T"] for r in rows]
    gt0 = [r["Brep_gt0_T"] for r in rows]
    gt001 = [r["Brep_gt001_T"] for r in rows]

    x = np.arange(len(labels), dtype=float)
    width = 0.25

    fig, ax = plt.subplots(figsize=(11.0, 6.2))

    ax.bar(x - width, geometry, width, label="Geometry selected; no threshold")
    ax.bar(x, gt0, width, label="|B| > 0")
    ax.bar(x + width, gt001, width, label="|B| > 0.01 T")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Representative |B| [T]")
    ax.set_title(
        "Threshold Sensitivity Check\n"
        "Primary result is geometry-selected; threshold versions are diagnostics"
    )
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_central_vs_full(path, full_sections, central_sections):
    labels = [s["section"] for s in full_sections]
    full_values = [s["Brep_geometry_T"] for s in full_sections]
    central_values = [s["Brep_geometry_T"] for s in central_sections]

    x = np.arange(len(labels), dtype=float)
    width = 0.34

    fig, ax = plt.subplots(figsize=(10.0, 6.0))

    ax.bar(x - width / 2, full_values, width, label="Full steel area")
    ax.bar(x + width / 2, central_values, width, label="Central/fiducial region")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Thickness-weighted representative |B| [T]")
    ax.set_title(
        "Full-Steel vs Central/Fiducial Representative Field\n"
        "Same plate-thickness weighting; no B threshold"
    )
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()

    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summary_table_lines(label, sections, full, momentum_gev):
    lines = []

    lines.append(label)
    lines.append("-" * len(label))

    for sec in sections:
        line = (
            f"{sec['section']:<13} "
            f"B_rep={sec['Brep_geometry_T']:.4f} T   "
            f"<B_perp>={sec['Bperp_rep_geometry_T']:.4f} T   "
            f"Integral(B_perp dl)={sec['integral_Bperp_Tm']:.4f} T*m"
        )

        if momentum_gev is not None:
            line += (
                f"   theta~{sec['approx_bend_angle_deg']:.3f} deg "
                f"at p={momentum_gev:g} GeV/c"
            )

        lines.append(line)

    lines.append("")
    lines.append(
        f"Full TMS     "
        f"B_rep={full['Brep_geometry_T']:.4f} T   "
        f"<B_perp>={full['Bperp_rep_geometry_T']:.4f} T   "
        f"Integral(B_perp dl)={full['integral_Bperp_Tm']:.4f} T*m"
    )

    if momentum_gev is not None:
        lines.append(
            f"Approx. straight-track bend at p={momentum_gev:g} GeV/c: "
            f"{full['approx_bend_angle_deg']:.3f} deg"
        )

    return lines


def write_summary(
    path,
    args,
    full_sections,
    full_summary,
    central_sections,
    central_summary,
    do_central,
    stream_stats,
):
    lines = []

    lines.append("GRAND TMS REPRESENTATIVE-FIELD ANALYSIS")
    lines.append("=" * 52)
    lines.append("")

    lines.append("PHYSICAL PLATE SECTIONS")
    lines.append("  Thin:         34 x 15 mm, starts 0 ... 2145 mm")
    lines.append("  Thick:        22 x 40 mm, starts 2210 ... 4100 mm")
    lines.append("  Double-thick: 24 x 80 mm, starts 4190 ... 7295 mm")
    lines.append("  Total steel thickness = 3310 mm")
    lines.append("")

    lines.append("PRIMARY DEFINITION")
    lines.append(
        "  Geometry decides which cells belong to steel; no B threshold."
    )
    lines.append(
        "  One mean |B| is calculated for each physical plate."
    )
    lines.append(
        "  Section/TMS representative values are weighted by real plate thickness."
    )
    lines.append("  No smoothing is used.")
    lines.append("")

    lines.append("FORMULA")
    lines.append(
        "  B_rep = sum_i( <|B|>_i * t_i ) / sum_i(t_i)"
    )
    lines.append("")

    lines.append("BENDING-POWER PROXY")
    lines.append(
        "  For beam-parallel +z reference: B_perp = sqrt(Bx^2 + By^2)"
    )
    lines.append(
        "  Integral(B_perp dl) ~= sum_i( <B_perp>_i * t_i )"
    )
    lines.append(
        "  This is a section-level bending-strength proxy, not yet a true "
        "track-by-track integral."
    )
    lines.append("")

    for line in summary_table_lines(
        "FULL STEEL AREA",
        full_sections,
        full_summary,
        args.momentum_gev,
    ):
        lines.append(line)

    lines.append("")
    lines.append("THRESHOLD DIAGNOSTICS")
    for sec in full_sections:
        lines.append(
            f"  {sec['section']:<13} "
            f"geometry={sec['Brep_geometry_T']:.4f} T, "
            f">0={sec['Brep_gt0_T']:.4f} T, "
            f">0.01={sec['Brep_gt001_T']:.4f} T, "
            f"mean zero fraction={100*sec['mean_zero_fraction']:.2f}%"
        )

    if do_central:
        lines.append("")
        lines.append(
            "CENTRAL/FIDUCIAL REGION "
            f"(|x|={args.central_x_min:g}...{args.central_x_max:g} mm, "
            f"|y|={args.central_y_min:g}...{args.central_y_max:g} mm)"
        )

        for line in summary_table_lines(
            "CENTRAL/FIDUCIAL",
            central_sections,
            central_summary,
            args.momentum_gev,
        ):
            lines.append(line)

    else:
        lines.append("")
        lines.append("CENTRAL/FIDUCIAL ANALYSIS")
        lines.append(
            "  Not run because no expert-defined central x-y bounds were supplied."
        )
        lines.append(
            "  Use --central-x-min, --central-x-max, --central-y-max "
            "(and optionally --central-y-min) when those bounds are agreed."
        )

    lines.append("")
    lines.append("STREAMING DIAGNOSTICS")
    lines.append(f"  Mapper data rows: {stream_stats['data_rows']:,}")
    lines.append(f"  XY steel-mask rows: {stream_stats['steel_rows']:,}")
    lines.append(f"  Rows inside physical plate z spans: {stream_stats['plate_rows']:,}")
    lines.append(f"  Rows in z gaps after XY mask: {stream_stats['gap_rows']:,}")

    if do_central:
        lines.append(
            f"  Central/fiducial plate rows: {stream_stats['central_rows']:,}"
        )

    lines.append("")
    lines.append("IMPORTANT CAVEATS")
    lines.append(
        "  * The transverse steel mask is still first-pass until exact "
        "GDML/CAD shoulder/opening dimensions are confirmed."
    )
    lines.append(
        "  * The current field dataset uses 135 mm pitch in the 80 mm section; "
        "the simulation geometry has been reported as 130 mm."
    )
    lines.append(
        "  * A true muon bend must integrate B_perp relative to the actual "
        "track direction and account for the evolving trajectory/momentum."
    )
    lines.append(
        "  * Do not tune the central region or B threshold to force the result "
        "toward an expected ~1 T."
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser():
    p = argparse.ArgumentParser(
        description="Grand corrected TMS representative-field analysis."
    )

    p.add_argument(
        "mapper",
        nargs="?",
        default="Mapper.txt",
    )

    p.add_argument(
        "-o",
        "--output-dir",
        default="plots/TMS_representative_B_grand",
    )

    # Diagnostic definition of exact/near zero.
    p.add_argument(
        "--zero-tolerance",
        type=float,
        default=1e-12,
    )

    # Optional muon momentum for a simple bend-angle estimate.
    p.add_argument(
        "--momentum-gev",
        type=float,
        default=None,
        help=(
            "Optional momentum [GeV/c]. If supplied, estimate bend angle from "
            "0.2998 * integral(B_perp dl) / p. This is only a straight-track proxy."
        ),
    )

    # First-pass transverse steel geometry.
    p.add_argument("--x-min", type=float, default=10.0)
    p.add_argument("--x-max", type=float, default=3730.0)
    p.add_argument("--y-max", type=float, default=2350.0)
    p.add_argument("--shoulder-y", type=float, default=1780.0)
    p.add_argument("--shoulder-x-min", type=float, default=250.0)
    p.add_argument("--shoulder-x-max", type=float, default=3480.0)
    p.add_argument("--opening-x-min", type=float, default=1620.0)
    p.add_argument("--opening-x-max", type=float, default=2120.0)
    p.add_argument("--opening-y-min", type=float, default=1780.0)
    p.add_argument("--opening-y-max", type=float, default=2000.0)

    # Optional expert-defined central/fiducial region.
    p.add_argument("--central-x-min", type=float, default=None)
    p.add_argument("--central-x-max", type=float, default=None)
    p.add_argument("--central-y-min", type=float, default=None)
    p.add_argument("--central-y-max", type=float, default=None)

    return p


def main():
    args = build_parser().parse_args()

    if args.momentum_gev is not None and args.momentum_gev <= 0:
        raise ValueError("--momentum-gev must be > 0.")

    mapper_path = Path(args.mapper)
    outdir = Path(args.output_dir)

    if not mapper_path.exists():
        raise FileNotFoundError(mapper_path)

    outdir.mkdir(parents=True, exist_ok=True)

    shape, header = read_mapper_header(mapper_path)
    ox, oy, oz, dx, dy, dz = header

    plates = make_plates(local_z0_mm=oz)
    verify_plate_definition(plates)

    do_central = validate_central_region(args)

    print("TMS plate definition")
    print("--------------------")
    print("Thin:         34 x 15 mm, starts 0 ... 2145 mm")
    print("Thick:        22 x 40 mm, starts 2210 ... 4100 mm")
    print("Double-thick: 24 x 80 mm, starts 4190 ... 7295 mm")
    print("Total steel thickness = 3310 mm")
    print()
    print(f"Mapper local source origin z0 = {oz:g} mm")
    print()
    print("Primary analysis:")
    print("  geometry-selected steel cells")
    print("  no B threshold")
    print("  no smoothing")
    print("  per-plate mean -> plate-thickness weighting")

    if do_central:
        print()
        print(
            "Central/fiducial analysis enabled:"
        )
        print(
            f"  |x| = {args.central_x_min:g} ... "
            f"{args.central_x_max:g} mm"
        )
        print(
            f"  |y| = {args.central_y_min:g} ... "
            f"{args.central_y_max:g} mm"
        )

    print()
    print("Reading Mapper.txt in one pass...")

    stream_stats = stream_mapper(
        mapper_path,
        plates,
        shape,
        header,
        args,
        do_central,
    )

    full_plate_rows = build_plate_rows(
        plates,
        stream_stats["full"],
    )

    full_sections = summarize_sections(full_plate_rows)
    full_summary = summarize_full(full_plate_rows)

    for sec in full_sections:
        add_bend_angle(sec, args.momentum_gev)

    add_bend_angle(full_summary, args.momentum_gev)

    full_plate_csv = outdir / "full_plate_results.csv"
    full_section_csv = outdir / "full_section_results.csv"

    write_csv(full_plate_csv, full_plate_rows)
    write_csv(
        full_section_csv,
        full_sections + [full_summary],
    )

    plot_plate_profile(
        outdir / "01_full_plate_profile.png",
        full_plate_rows,
        full_sections,
        full_summary,
        "TMS Representative Magnetic Field by Physical Plate Thickness",
    )

    plot_bending_power(
        outdir / "02_full_bending_power_by_section.png",
        full_sections,
        "Full Steel Area",
    )

    plot_threshold_sensitivity(
        outdir / "03_threshold_sensitivity.png",
        full_sections,
        full_summary,
    )

    central_plate_rows = None
    central_sections = None
    central_summary = None

    if do_central:
        central_plate_rows = build_plate_rows(
            plates,
            stream_stats["central"],
        )

        central_sections = summarize_sections(
            central_plate_rows,
        )

        central_summary = summarize_full(
            central_plate_rows,
        )

        for sec in central_sections:
            add_bend_angle(sec, args.momentum_gev)

        add_bend_angle(
            central_summary,
            args.momentum_gev,
        )

        write_csv(
            outdir / "central_plate_results.csv",
            central_plate_rows,
        )

        write_csv(
            outdir / "central_section_results.csv",
            central_sections + [central_summary],
        )

        plot_central_vs_full(
            outdir / "04_central_vs_full_sections.png",
            full_sections,
            central_sections,
        )

        plot_plate_profile(
            outdir / "05_central_plate_profile.png",
            central_plate_rows,
            central_sections,
            central_summary,
            "Central/Fiducial TMS Representative Field",
        )

        plot_bending_power(
            outdir / "06_central_bending_power_by_section.png",
            central_sections,
            "Central/Fiducial Region",
        )

    write_summary(
        outdir / "grand_summary.txt",
        args,
        full_sections,
        full_summary,
        central_sections,
        central_summary,
        do_central,
        stream_stats,
    )

    print()
    print("FULL-STEEL RESULTS")
    print("------------------")

    for sec in full_sections:
        print(
            f"{sec['section']:<13} "
            f"B_rep={sec['Brep_geometry_T']:.4f} T   "
            f"Integral(B_perp dl)={sec['integral_Bperp_Tm']:.4f} T*m"
        )

    print()
    print(
        f"Full TMS B_rep = {full_summary['Brep_geometry_T']:.4f} T"
    )
    print(
        f"Full TMS Integral(B_perp dl) = "
        f"{full_summary['integral_Bperp_Tm']:.4f} T*m"
    )

    if args.momentum_gev is not None:
        print(
            f"Straight-track bend proxy at p={args.momentum_gev:g} GeV/c: "
            f"{full_summary['approx_bend_angle_deg']:.3f} deg"
        )

    print()
    print("Threshold sensitivity (diagnostic only):")
    for sec in full_sections:
        print(
            f"  {sec['section']:<13} "
            f"geometry={sec['Brep_geometry_T']:.4f} T   "
            f">0={sec['Brep_gt0_T']:.4f} T   "
            f">0.01={sec['Brep_gt001_T']:.4f} T"
        )

    if do_central:
        print()
        print("CENTRAL/FIDUCIAL RESULTS")
        print("------------------------")

        for sec in central_sections:
            print(
                f"{sec['section']:<13} "
                f"B_rep={sec['Brep_geometry_T']:.4f} T   "
                f"Integral(B_perp dl)={sec['integral_Bperp_Tm']:.4f} T*m"
            )

        print()
        print(
            f"Central TMS B_rep = "
            f"{central_summary['Brep_geometry_T']:.4f} T"
        )

    print()
    print("Saved to:")
    print(f"  {outdir}")
    print()
    print("Most important files:")
    print(f"  {outdir / '01_full_plate_profile.png'}")
    print(f"  {outdir / '02_full_bending_power_by_section.png'}")
    print(f"  {outdir / '03_threshold_sensitivity.png'}")
    print(f"  {outdir / 'full_section_results.csv'}")
    print(f"  {outdir / 'grand_summary.txt'}")

    if not do_central:
        print()
        print(
            "Central/fiducial comparison was NOT run because no x-y bounds "
            "were supplied."
        )
        print(
            "Once Marco/expert provides the central region, rerun with:"
        )
        print(
            "  --central-x-min XMIN --central-x-max XMAX "
            "--central-y-max YMAX"
        )


if __name__ == "__main__":
    main()
