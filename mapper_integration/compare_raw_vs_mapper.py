#!/usr/bin/env python3
"""Validate Mapper.txt against the original .fld magnetic-field values.

This script intentionally keeps mapper validation separate from the later
representative-field analysis.

Validation rules
----------------
- No |B| > 0.01 T cut is applied here.
- The original .fld values are treated as the source magnetic-field data.
- Mapper.txt is evaluated back at the original source-point coordinates using
  plate-aware trilinear interpolation on the uniform mapper grid.
- The same four visual comparisons are produced at one selected z layer:

    1. quarter geometry: quiver vs quiver
    2. quarter geometry: |B| heatmap vs |B| heatmap
    3. full geometry: quiver vs quiver
    4. full geometry: |B| heatmap vs |B| heatmap

- Numerical comparison tables are also produced for Bx, By, Bz, |B| and for
  the working longitudinal TMS regions.

The raw heatmaps are formed by binning the original irregular .fld samples
onto the Mapper.txt XY cell centres for display only. Empty raw cells stay
blank; they are not treated as B=0.

The numerical validation does NOT use those binned heatmaps. It compares each
original raw source point directly with Mapper.txt evaluated at that same
physical coordinate.

Quiver-plot convention
----------------------
Arrow lengths are normalized to emphasize field direction. Arrow colour carries
the magnetic-field magnitude |B|. Therefore a visible arrow near a boundary
does not imply that its field magnitude is as large as an arrow in the steel.

Masking convention
------------------
The current mapper uses a source-support radius to suppress unsupported
extrapolation. This is NOT yet a true steel-geometry mask. A future
geometry-aware mask should use the exact latest transverse steel outline,
cut-outs, and gaps. The current validation therefore tests the mapper exactly
as it is presently used in simulation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from tms_mapper.core import determine_z_placement, discover_plate_groups, read_fld
from tms_mapper.edep import iter_mapper_rows, read_mapper_header


DEFAULT_REGION1_END_MM = -500.0
DEFAULT_REGION2_END_MM = 2000.0
DEFAULT_REGION3_END_MM = 2900.0


@dataclass
class ValidationArrays:
    points: np.ndarray
    raw_vectors: np.ndarray
    mapper_vectors: np.ndarray


def load_mapper_vector_cube(
    mapper_path: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, object]:
    """Load Bx, By, Bz into an (NX, NY, NZ, 3) cube."""
    header = read_mapper_header(mapper_path)
    if header.shape is None:
        raise ValueError("Mapper.txt needs '# shape NX NY NZ' metadata.")

    nx, ny, nz = header.shape
    expected = nx * ny * nz
    vectors = np.empty((expected, 3), dtype=np.float32)

    count = 0
    for row in iter_mapper_rows(mapper_path):
        if count >= expected:
            raise ValueError("Mapper has more rows than declared by '# shape'.")
        vectors[count] = row[3:6]
        count += 1

    if count != expected:
        raise ValueError(
            f"Mapper row count {count:,} != declared grid size {expected:,}."
        )

    x_mm = header.offset_x_mm + np.arange(nx, dtype=float) * header.dx_mm
    y_mm = header.offset_y_mm + np.arange(ny, dtype=float) * header.dy_mm
    z_mm = header.offset_z_mm + np.arange(nz, dtype=float) * header.dz_mm

    return vectors.reshape(nx, ny, nz, 3), x_mm, y_mm, z_mm, header


def aligned_group_data(
    group,
    *,
    source_length_unit: str,
    z_placement: str,
    base_z_mm: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Load Plate1+Plate2 and apply the same z placement used by the mapper."""
    point_parts: List[np.ndarray] = []
    field_parts: List[np.ndarray] = []
    expected_start = base_z_mm + group.z_offset_mm

    for part in (1, 2):
        points, fields = read_fld(
            group.maps[part].path,
            source_length_unit=source_length_unit,
        )
        if z_placement == "zcoord":
            shift = expected_start - float(np.min(points[:, 2]))
            points = points.copy()
            points[:, 2] += shift

        point_parts.append(points)
        field_parts.append(fields)

    return np.concatenate(point_parts), np.concatenate(field_parts)


def select_example_layer(
    groups,
    *,
    example_z_mm: float,
    source_length_unit: str,
    z_placement: str,
    base_z_mm: float,
):
    """Choose the physical plate and raw source z-layer nearest requested z."""
    best = None

    for group in groups:
        points, fields = aligned_group_data(
            group,
            source_length_unit=source_length_unit,
            z_placement=z_placement,
            base_z_mm=base_z_mm,
        )
        zmin = float(points[:, 2].min())
        zmax = float(points[:, 2].max())

        if zmin <= example_z_mm <= zmax:
            distance = 0.0
        else:
            distance = min(abs(example_z_mm - zmin), abs(example_z_mm - zmax))

        if best is None or distance < best[0]:
            best = (distance, group, points, fields)

    if best is None:
        raise ValueError("No source plate groups were found.")

    _, group, points, fields = best
    unique_z = np.unique(np.round(points[:, 2], 6))
    raw_z = float(unique_z[np.argmin(np.abs(unique_z - example_z_mm))])

    mask = np.isclose(points[:, 2], raw_z, atol=1e-5)
    return group, points[mask], fields[mask], raw_z


def reflect_quarter_to_full(
    points: np.ndarray,
    fields: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Reflect +X,+Y source data using the mapper axial-vector B symmetry."""
    p_parts: List[np.ndarray] = []
    b_parts: List[np.ndarray] = []

    for sx, sy in ((1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)):
        p = points.copy()
        p[:, 0] *= sx
        p[:, 1] *= sy

        b = fields.copy()
        b[:, 0] *= sy
        b[:, 1] *= sx
        b[:, 2] *= sx * sy

        p_parts.append(p)
        b_parts.append(b)

    return np.concatenate(p_parts), np.concatenate(b_parts)


def _check_points_in_xy_bounds(
    points: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
) -> None:
    tol_x = 1e-8 + 0.5 * abs(float(x_axis[1] - x_axis[0]))
    tol_y = 1e-8 + 0.5 * abs(float(y_axis[1] - y_axis[0]))

    if (
        np.any(points[:, 0] < x_axis[0] - tol_x)
        or np.any(points[:, 0] > x_axis[-1] + tol_x)
        or np.any(points[:, 1] < y_axis[0] - tol_y)
        or np.any(points[:, 1] > y_axis[-1] + tol_y)
    ):
        raise ValueError("Some source points lie outside Mapper.txt XY bounds.")


def _lower_index_and_fraction(
    values: np.ndarray,
    axis: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return lower grid index and interpolation fraction for a uniform axis."""
    if len(axis) < 2:
        raise ValueError("Axis needs at least two grid points.")

    spacing = float(axis[1] - axis[0])
    pos = (values - axis[0]) / spacing
    lower = np.floor(pos).astype(int)
    lower = np.clip(lower, 0, len(axis) - 2)

    frac = (values - axis[lower]) / spacing
    frac = np.clip(frac, 0.0, 1.0)
    return lower, frac


def _bilinear_xy_at_indices(
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    points: np.ndarray,
    iz: np.ndarray,
) -> np.ndarray:
    """Bilinearly interpolate XY for per-point mapper z indices."""
    ix, fx = _lower_index_and_fraction(points[:, 0], x_axis)
    iy, fy = _lower_index_and_fraction(points[:, 1], y_axis)

    v00 = cube[ix, iy, iz]
    v10 = cube[ix + 1, iy, iz]
    v01 = cube[ix, iy + 1, iz]
    v11 = cube[ix + 1, iy + 1, iz]

    fx = fx[:, None]
    fy = fy[:, None]

    return (
        (1.0 - fx) * (1.0 - fy) * v00
        + fx * (1.0 - fy) * v10
        + (1.0 - fx) * fy * v01
        + fx * fy * v11
    ).astype(np.float64)


def sample_mapper_at_source_points(
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    points: np.ndarray,
    *,
    plate_zmin: float,
    plate_zmax: float,
) -> np.ndarray:
    """Evaluate Mapper.txt at original source coordinates.

    X and Y use bilinear interpolation.

    Z uses linear interpolation only between mapper z planes that lie inside
    the same physical source plate slab. This prevents a source point at the
    edge of a plate from being numerically mixed with a zero-filled detector
    gap outside that plate.

    If a source z lies beyond the first/last mapper plane contained in the
    source plate slab, the nearest in-slab mapper plane is used at that edge.
    """
    _check_points_in_xy_bounds(points, x_axis, y_axis)

    valid_z_idx = np.where(
        (z_axis >= plate_zmin - 1e-6) & (z_axis <= plate_zmax + 1e-6)
    )[0]

    if valid_z_idx.size == 0:
        raise ValueError(
            f"No mapper z plane lies inside source plate slab "
            f"[{plate_zmin:g}, {plate_zmax:g}] mm."
        )

    if valid_z_idx.size == 1:
        iz = np.full(len(points), valid_z_idx[0], dtype=int)
        return _bilinear_xy_at_indices(cube, x_axis, y_axis, points, iz)

    valid_z = z_axis[valid_z_idx]
    source_z = points[:, 2]

    insertion = np.searchsorted(valid_z, source_z, side="right")
    upper_local = np.clip(insertion, 0, len(valid_z) - 1)
    lower_local = np.clip(insertion - 1, 0, len(valid_z) - 1)

    below = source_z <= valid_z[0]
    above = source_z >= valid_z[-1]
    lower_local[below] = 0
    upper_local[below] = 0
    lower_local[above] = len(valid_z) - 1
    upper_local[above] = len(valid_z) - 1

    iz0 = valid_z_idx[lower_local]
    iz1 = valid_z_idx[upper_local]

    b0 = _bilinear_xy_at_indices(cube, x_axis, y_axis, points, iz0)
    b1 = _bilinear_xy_at_indices(cube, x_axis, y_axis, points, iz1)

    z0 = z_axis[iz0]
    z1 = z_axis[iz1]

    denom = z1 - z0
    fz = np.zeros(len(points), dtype=float)
    interp = np.abs(denom) > 1e-12
    fz[interp] = (source_z[interp] - z0[interp]) / denom[interp]
    fz = np.clip(fz, 0.0, 1.0)[:, None]

    return (1.0 - fz) * b0 + fz * b1


def error_summary(raw: np.ndarray, mapped: np.ndarray) -> Dict[str, float]:
    """Return expert-facing agreement metrics for one scalar quantity."""
    raw = np.asarray(raw, dtype=float)
    mapped = np.asarray(mapped, dtype=float)

    delta = mapped - raw
    abs_delta = np.abs(delta)

    if raw.size > 1 and np.std(raw) > 0.0 and np.std(mapped) > 0.0:
        corr = float(np.corrcoef(raw, mapped)[0, 1])
    else:
        corr = float("nan")

    return {
        "n_points": int(raw.size),
        "raw_mean_T": float(np.mean(raw)),
        "mapper_mean_T": float(np.mean(mapped)),
        "bias_T": float(np.mean(delta)),
        "mae_T": float(np.mean(abs_delta)),
        "rmse_T": float(np.sqrt(np.mean(delta * delta))),
        "p95_abs_error_T": float(np.percentile(abs_delta, 95)),
        "max_abs_error_T": float(np.max(abs_delta)),
        "correlation": corr,
    }


def _normalised_xy_arrows(fields: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    u = np.asarray(fields[:, 0], dtype=float)
    v = np.asarray(fields[:, 1], dtype=float)
    norm = np.hypot(u, v)

    out_u = np.zeros_like(u)
    out_v = np.zeros_like(v)

    good = norm > 1e-12
    out_u[good] = u[good] / norm[good]
    out_v[good] = v[good] / norm[good]
    return out_u, out_v


def _subsample_indices(n: int, max_points: int) -> np.ndarray:
    if n <= max_points:
        return np.arange(n, dtype=int)
    return np.linspace(0, n - 1, max_points, dtype=int)


def _mapper_slice_points(
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    iz: int,
    *,
    quarter: bool,
    step: int,
) -> Tuple[np.ndarray, np.ndarray]:
    if quarter:
        ix = np.where(x_axis >= -1e-9)[0][::step]
        iy = np.where(y_axis >= -1e-9)[0][::step]
    else:
        ix = np.arange(0, len(x_axis), step, dtype=int)
        iy = np.arange(0, len(y_axis), step, dtype=int)

    xx, yy = np.meshgrid(x_axis[ix], y_axis[iy], indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), np.full(xx.size, np.nan)))
    fields = cube[np.ix_(ix, iy, [iz])][:, :, 0, :].reshape(-1, 3)
    return points, fields


def _quiver_panel(
    ax,
    points: np.ndarray,
    fields: np.ndarray,
    *,
    max_points: int,
    vmin: float,
    vmax: float,
    title: str,
):
    idx = _subsample_indices(len(points), max_points)
    p = points[idx]
    b = fields[idx]

    mag = np.linalg.norm(b, axis=1)
    u, v = _normalised_xy_arrows(b)

    q = ax.quiver(
        p[:, 0],
        p[:, 1],
        u,
        v,
        mag,
        angles="xy",
        scale_units="xy",
        scale=0.0065,
        width=0.0025,
        pivot="mid",
        clim=(vmin, vmax),
    )

    ax.set_title(title)
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.12)
    return q


def _set_quarter_limits(axes, raw_points: np.ndarray, x_axis, y_axis) -> None:
    x_max = max(
        float(np.max(raw_points[:, 0])),
        float(np.max(x_axis[x_axis >= 0.0])),
    )
    y_max = max(
        float(np.max(raw_points[:, 1])),
        float(np.max(y_axis[y_axis >= 0.0])),
    )

    for ax in axes:
        ax.set_xlim(0.0, x_max)
        ax.set_ylim(0.0, y_max)


def _set_full_limits(axes, raw_points: np.ndarray, x_axis, y_axis) -> None:
    x_max = max(
        float(np.max(np.abs(raw_points[:, 0]))),
        float(np.max(np.abs(x_axis))),
    )
    y_max = max(
        float(np.max(np.abs(raw_points[:, 1]))),
        float(np.max(np.abs(y_axis))),
    )

    for ax in axes:
        ax.set_xlim(-x_max, x_max)
        ax.set_ylim(-y_max, y_max)


def make_quarter_quiver_comparison(
    path: Path,
    raw_points: np.ndarray,
    raw_fields: np.ndarray,
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    raw_z: float,
    *,
    quiver_step: int,
    raw_quiver_max: int,
    vmin: float,
    vmax: float,
) -> float:
    iz = int(np.argmin(np.abs(z_axis - raw_z)))
    mapper_points, mapper_fields = _mapper_slice_points(
        cube, x_axis, y_axis, iz, quarter=True, step=quiver_step
    )

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2), constrained_layout=True)

    q = _quiver_panel(
        axes[0],
        raw_points,
        raw_fields,
        max_points=raw_quiver_max,
        vmin=vmin,
        vmax=vmax,
        title=f"Original .fld vectors\nz ≈ {raw_z:g} mm",
    )

    _quiver_panel(
        axes[1],
        mapper_points,
        mapper_fields,
        max_points=len(mapper_points),
        vmin=vmin,
        vmax=vmax,
        title=f"Mapper.txt quarter vectors\nz = {z_axis[iz]:g} mm",
    )

    _set_quarter_limits(axes, raw_points, x_axis, y_axis)

    fig.suptitle(
        "Quarter geometry: original .fld vs Mapper.txt — field direction\n"
        "Arrow direction = in-plane field direction; colour = |B|",
        fontsize=14,
    )

    cbar = fig.colorbar(q, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("|B| [T]")

    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)

    return float(z_axis[iz])


def make_full_quiver_comparison(
    path: Path,
    raw_quarter_points: np.ndarray,
    raw_quarter_fields: np.ndarray,
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    raw_z: float,
    *,
    quiver_step: int,
    raw_quiver_max: int,
    vmin: float,
    vmax: float,
) -> float:
    raw_points, raw_fields = reflect_quarter_to_full(
        raw_quarter_points,
        raw_quarter_fields,
    )

    iz = int(np.argmin(np.abs(z_axis - raw_z)))
    mapper_points, mapper_fields = _mapper_slice_points(
        cube, x_axis, y_axis, iz, quarter=False, step=quiver_step
    )

    fig, axes = plt.subplots(1, 2, figsize=(15, 7.0), constrained_layout=True)

    q = _quiver_panel(
        axes[0],
        raw_points,
        raw_fields,
        max_points=raw_quiver_max,
        vmin=vmin,
        vmax=vmax,
        title=f"Original .fld quarter reflected to full geometry\nz ≈ {raw_z:g} mm",
    )

    _quiver_panel(
        axes[1],
        mapper_points,
        mapper_fields,
        max_points=len(mapper_points),
        vmin=vmin,
        vmax=vmax,
        title=f"Full Mapper.txt\nz = {z_axis[iz]:g} mm",
    )

    _set_full_limits(axes, raw_points, x_axis, y_axis)

    fig.suptitle(
        "Full geometry: reflected .fld vs Mapper.txt — field direction\n"
        "Arrow direction = in-plane field direction; colour = |B|",
        fontsize=14,
    )

    cbar = fig.colorbar(q, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("|B| [T]")

    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)

    return float(z_axis[iz])


def _raw_binned_heatmap(
    points: np.ndarray,
    fields: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
) -> np.ndarray:
    """Bin irregular raw samples onto mapper XY cell centres for display only."""
    dx = float(x_axis[1] - x_axis[0])
    dy = float(y_axis[1] - y_axis[0])

    ix = np.rint((points[:, 0] - x_axis[0]) / dx).astype(int)
    iy = np.rint((points[:, 1] - y_axis[0]) / dy).astype(int)

    valid = (
        (ix >= 0)
        & (ix < len(x_axis))
        & (iy >= 0)
        & (iy < len(y_axis))
    )

    ix = ix[valid]
    iy = iy[valid]
    mag = np.linalg.norm(fields[valid], axis=1)

    sums = np.zeros((len(x_axis), len(y_axis)), dtype=float)
    counts = np.zeros((len(x_axis), len(y_axis)), dtype=int)

    np.add.at(sums, (ix, iy), mag)
    np.add.at(counts, (ix, iy), 1)

    result = np.full_like(sums, np.nan, dtype=float)
    populated = counts > 0
    result[populated] = sums[populated] / counts[populated]
    return result


def _heatmap_pair(
    path: Path,
    raw_points: np.ndarray,
    raw_fields: np.ndarray,
    mapper_mag: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    *,
    quarter: bool,
    raw_z: float,
    mapper_z: float,
    vmin: float,
    vmax: float,
    title: str,
) -> None:
    raw_grid = _raw_binned_heatmap(raw_points, raw_fields, x_axis, y_axis)

    if quarter:
        ix = np.where(x_axis >= -1e-9)[0]
        iy = np.where(y_axis >= -1e-9)[0]
    else:
        ix = np.arange(len(x_axis), dtype=int)
        iy = np.arange(len(y_axis), dtype=int)

    raw_view = raw_grid[np.ix_(ix, iy)]
    mapper_view = mapper_mag[np.ix_(ix, iy)]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4), constrained_layout=True)

    axes[0].pcolormesh(
        x_axis[ix],
        y_axis[iy],
        np.ma.masked_invalid(raw_view).T,
        shading="nearest",
        vmin=vmin,
        vmax=vmax,
    )

    mapper_mesh = axes[1].pcolormesh(
        x_axis[ix],
        y_axis[iy],
        mapper_view.T,
        shading="nearest",
        vmin=vmin,
        vmax=vmax,
    )

    for ax in axes:
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("X [mm]")
        ax.set_ylabel("Y [mm]")
        ax.grid(True, alpha=0.12)

    axes[0].set_title(
        f"Original .fld |B| heatmap\n"
        f"z ≈ {raw_z:g} mm; blank = no raw sample in that cell"
    )
    axes[1].set_title(f"Mapper.txt |B| heatmap\nz = {mapper_z:g} mm")

    if quarter:
        _set_quarter_limits(axes, raw_points, x_axis, y_axis)
    else:
        _set_full_limits(axes, raw_points, x_axis, y_axis)

    fig.suptitle(title, fontsize=15)

    cbar = fig.colorbar(mapper_mesh, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("|B| [T]")

    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def make_quarter_heatmap_comparison(
    path: Path,
    raw_points: np.ndarray,
    raw_fields: np.ndarray,
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    raw_z: float,
    *,
    vmin: float,
    vmax: float,
) -> float:
    iz = int(np.argmin(np.abs(z_axis - raw_z)))
    mapper_mag = np.linalg.norm(cube[:, :, iz, :], axis=2)

    _heatmap_pair(
        path,
        raw_points,
        raw_fields,
        mapper_mag,
        x_axis,
        y_axis,
        quarter=True,
        raw_z=raw_z,
        mapper_z=float(z_axis[iz]),
        vmin=vmin,
        vmax=vmax,
        title="Quarter geometry: original .fld vs Mapper.txt — |B|",
    )

    return float(z_axis[iz])


def make_full_heatmap_comparison(
    path: Path,
    raw_quarter_points: np.ndarray,
    raw_quarter_fields: np.ndarray,
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    raw_z: float,
    *,
    vmin: float,
    vmax: float,
) -> float:
    raw_points, raw_fields = reflect_quarter_to_full(
        raw_quarter_points,
        raw_quarter_fields,
    )

    iz = int(np.argmin(np.abs(z_axis - raw_z)))
    mapper_mag = np.linalg.norm(cube[:, :, iz, :], axis=2)

    _heatmap_pair(
        path,
        raw_points,
        raw_fields,
        mapper_mag,
        x_axis,
        y_axis,
        quarter=False,
        raw_z=raw_z,
        mapper_z=float(z_axis[iz]),
        vmin=vmin,
        vmax=vmax,
        title="Full geometry: reflected .fld vs Mapper.txt — |B|",
    )

    return float(z_axis[iz])


def collect_numerical_validation(
    groups,
    cube: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    *,
    source_length_unit: str,
    z_placement: str,
    base_z_mm: float,
) -> ValidationArrays:
    """Compare all original quarter source points to Mapper.txt at same XYZ."""
    point_parts: List[np.ndarray] = []
    raw_parts: List[np.ndarray] = []
    mapper_parts: List[np.ndarray] = []

    for index, group in enumerate(groups, 1):
        points, raw_vectors = aligned_group_data(
            group,
            source_length_unit=source_length_unit,
            z_placement=z_placement,
            base_z_mm=base_z_mm,
        )

        zmin = float(points[:, 2].min())
        zmax = float(points[:, 2].max())

        mapper_vectors = sample_mapper_at_source_points(
            cube,
            x_axis,
            y_axis,
            z_axis,
            points,
            plate_zmin=zmin,
            plate_zmax=zmax,
        )

        point_parts.append(points)
        raw_parts.append(raw_vectors)
        mapper_parts.append(mapper_vectors)

        if index % 10 == 0 or index == len(groups):
            print(f"  compared {index}/{len(groups)} physical plate positions")

    return ValidationArrays(
        points=np.concatenate(point_parts),
        raw_vectors=np.concatenate(raw_parts),
        mapper_vectors=np.concatenate(mapper_parts),
    )


def build_component_rows(data: ValidationArrays) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []

    raw_mag = np.linalg.norm(data.raw_vectors, axis=1)
    mapper_mag = np.linalg.norm(data.mapper_vectors, axis=1)

    quantities = [
        ("Bx", data.raw_vectors[:, 0], data.mapper_vectors[:, 0]),
        ("By", data.raw_vectors[:, 1], data.mapper_vectors[:, 1]),
        ("Bz", data.raw_vectors[:, 2], data.mapper_vectors[:, 2]),
        ("|B|", raw_mag, mapper_mag),
    ]

    for name, raw, mapped in quantities:
        row = {"quantity": name}
        row.update(error_summary(raw, mapped))
        rows.append(row)

    return rows


def _region_mask(
    z: np.ndarray,
    zlo: float,
    zhi: float,
    *,
    include_upper: bool = False,
) -> np.ndarray:
    if include_upper:
        return (z >= zlo) & (z <= zhi)
    return (z >= zlo) & (z < zhi)


def build_region_rows(
    data: ValidationArrays,
    *,
    region1_end_mm: float,
    region2_end_mm: float,
    region3_end_mm: float,
) -> List[Dict[str, float]]:
    raw_mag = np.linalg.norm(data.raw_vectors, axis=1)
    mapper_mag = np.linalg.norm(data.mapper_vectors, axis=1)
    z = data.points[:, 2]

    zmin = float(np.min(z))
    zmax = float(np.max(z))

    definitions = [
        ("Region 1", zmin, region1_end_mm, False),
        ("Region 2", region1_end_mm, region2_end_mm, False),
        ("Region 3", region2_end_mm, region3_end_mm, True),
        ("Whole TMS", zmin, zmax, True),
    ]

    rows: List[Dict[str, float]] = []

    for name, zlo, zhi, include_upper in definitions:
        mask = _region_mask(z, zlo, zhi, include_upper=include_upper)
        if not np.any(mask):
            continue

        stats = error_summary(raw_mag[mask], mapper_mag[mask])

        row = {
            "region": name,
            "z_start_mm": float(zlo),
            "z_end_mm": float(zhi),
            "n_points": int(np.sum(mask)),
            "raw_mean_B_T": stats["raw_mean_T"],
            "mapper_mean_B_T": stats["mapper_mean_T"],
            "difference_T": stats["bias_T"],
            "mae_T": stats["mae_T"],
            "rmse_T": stats["rmse_T"],
            "p95_abs_error_T": stats["p95_abs_error_T"],
            "max_abs_error_T": stats["max_abs_error_T"],
            "correlation": stats["correlation"],
        }
        rows.append(row)

    return rows


def write_component_csv(path: Path, rows: Sequence[Dict[str, float]]) -> None:
    fields = [
        "quantity",
        "n_points",
        "raw_mean_T",
        "mapper_mean_T",
        "bias_T",
        "mae_T",
        "rmse_T",
        "p95_abs_error_T",
        "max_abs_error_T",
        "correlation",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_region_csv(path: Path, rows: Sequence[Dict[str, float]]) -> None:
    fields = [
        "region",
        "z_start_mm",
        "z_end_mm",
        "n_points",
        "raw_mean_B_T",
        "mapper_mean_B_T",
        "difference_T",
        "mae_T",
        "rmse_T",
        "p95_abs_error_T",
        "max_abs_error_T",
        "correlation",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _fmt(value, digits: int = 5) -> str:
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"

    value = float(value)
    if np.isnan(value):
        return "n/a"
    return f"{value:.{digits}f}"


def make_component_table_png(
    path: Path,
    rows: Sequence[Dict[str, float]],
) -> None:
    columns = [
        ("Quantity", "quantity"),
        ("N", "n_points"),
        ("Raw mean\n[T]", "raw_mean_T"),
        ("Mapper mean\n[T]", "mapper_mean_T"),
        ("Bias\n[T]", "bias_T"),
        ("MAE\n[T]", "mae_T"),
        ("RMSE\n[T]", "rmse_T"),
        ("95% abs.\nerror [T]", "p95_abs_error_T"),
        ("Max abs.\nerror [T]", "max_abs_error_T"),
        ("Correlation", "correlation"),
    ]

    table_data = []
    for row in rows:
        formatted = []
        for _, key in columns:
            if key == "quantity":
                formatted.append(str(row[key]))
            elif key == "n_points":
                formatted.append(_fmt(row[key], 0))
            elif key == "correlation":
                formatted.append(_fmt(row[key], 6))
            else:
                formatted.append(_fmt(row[key], 6))
        table_data.append(formatted)

    fig, ax = plt.subplots(figsize=(17, 3.3))
    ax.axis("off")

    table = ax.table(
        cellText=table_data,
        colLabels=[label for label, _ in columns],
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.75)

    ax.set_title(
        "Raw .fld vs Mapper.txt numerical validation at matching source coordinates\n"
        "No magnetic-field threshold applied",
        pad=16,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def make_region_table_png(
    path: Path,
    rows: Sequence[Dict[str, float]],
) -> None:
    columns = [
        ("Region", "region"),
        ("Z range\n[mm]", "z_range"),
        ("N", "n_points"),
        ("Raw mean |B|\n[T]", "raw_mean_B_T"),
        ("Mapper mean |B|\n[T]", "mapper_mean_B_T"),
        ("Difference\n[T]", "difference_T"),
        ("MAE\n[T]", "mae_T"),
        ("RMSE\n[T]", "rmse_T"),
        ("95% abs.\nerror [T]", "p95_abs_error_T"),
        ("Correlation", "correlation"),
    ]

    table_data = []
    for row in rows:
        display = dict(row)
        display["z_range"] = f"{row['z_start_mm']:.0f} to {row['z_end_mm']:.0f}"

        formatted = []
        for _, key in columns:
            if key in {"region", "z_range"}:
                formatted.append(str(display[key]))
            elif key == "n_points":
                formatted.append(_fmt(display[key], 0))
            elif key == "correlation":
                formatted.append(_fmt(display[key], 6))
            else:
                formatted.append(_fmt(display[key], 6))
        table_data.append(formatted)

    fig, ax = plt.subplots(figsize=(17, 3.5))
    ax.axis("off")

    table = ax.table(
        cellText=table_data,
        colLabels=[label for label, _ in columns],
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.75)

    ax.set_title(
        "Raw .fld vs Mapper.txt |B| agreement by working TMS region\n"
        "No magnetic-field threshold applied",
        pad=16,
    )

    fig.tight_layout()
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def print_component_table(rows: Sequence[Dict[str, float]]) -> None:
    print("\nNumerical validation at original .fld source coordinates")
    print("No magnetic-field threshold applied.")
    print(
        f"{'Quantity':<9} {'N':>10} {'Raw mean':>11} {'Map mean':>11} "
        f"{'Bias':>11} {'MAE':>11} {'RMSE':>11} {'P95':>11} {'Corr':>10}"
    )
    print("-" * 109)

    for row in rows:
        print(
            f"{row['quantity']:<9} "
            f"{int(row['n_points']):>10,} "
            f"{row['raw_mean_T']:>11.6f} "
            f"{row['mapper_mean_T']:>11.6f} "
            f"{row['bias_T']:>+11.6f} "
            f"{row['mae_T']:>11.6f} "
            f"{row['rmse_T']:>11.6f} "
            f"{row['p95_abs_error_T']:>11.6f} "
            f"{row['correlation']:>10.6f}"
        )


def print_region_table(rows: Sequence[Dict[str, float]]) -> None:
    print("\n|B| agreement by working TMS region")
    print(
        f"{'Region':<10} {'Z range [mm]':>20} {'N':>10} "
        f"{'Raw mean':>11} {'Map mean':>11} {'Diff':>11} "
        f"{'MAE':>11} {'RMSE':>11}"
    )
    print("-" * 111)

    for row in rows:
        z_range = f"{row['z_start_mm']:.0f}..{row['z_end_mm']:.0f}"
        print(
            f"{row['region']:<10} "
            f"{z_range:>20} "
            f"{int(row['n_points']):>10,} "
            f"{row['raw_mean_B_T']:>11.6f} "
            f"{row['mapper_mean_B_T']:>11.6f} "
            f"{row['difference_T']:>+11.6f} "
            f"{row['mae_T']:>11.6f} "
            f"{row['rmse_T']:>11.6f}"
        )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate four raw-.fld-vs-Mapper.txt validation plots and "
            "matching-coordinate numerical validation tables."
        )
    )

    parser.add_argument("mapper", nargs="?", default="Mapper.txt")
    parser.add_argument("--field-dir", default="Field_maps")
    parser.add_argument(
        "--z",
        "--example-z",
        dest="example_z",
        type=float,
        default=-4000.0,
        help="Requested local-z slice for the four visual comparison plots.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="plots/raw_vs_mapper_validation",
    )
    parser.add_argument(
        "--source-length-unit",
        choices=("m", "mm"),
        default="m",
    )
    parser.add_argument(
        "--z-placement",
        choices=("auto", "source", "zcoord"),
        default="auto",
    )
    parser.add_argument("--quiver-step", type=int, default=2)
    parser.add_argument("--raw-quiver-max", type=int, default=1400)

    parser.add_argument(
        "--region1-end",
        type=float,
        default=DEFAULT_REGION1_END_MM,
    )
    parser.add_argument(
        "--region2-end",
        type=float,
        default=DEFAULT_REGION2_END_MM,
    )
    parser.add_argument(
        "--region3-end",
        type=float,
        default=DEFAULT_REGION3_END_MM,
    )

    return parser


def main() -> int:
    args = make_parser().parse_args()

    if args.quiver_step < 1:
        raise ValueError("--quiver-step must be >= 1.")
    if args.raw_quiver_max < 1:
        raise ValueError("--raw-quiver-max must be positive.")

    mapper_path = Path(args.mapper)
    field_dir = Path(args.field_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading Mapper.txt ...")
    cube, x_axis, y_axis, z_axis, header = load_mapper_vector_cube(mapper_path)

    print("Reading original .fld plate groups ...")
    groups = discover_plate_groups(field_dir)
    resolved_z_mode, base_z_mm = determine_z_placement(
        groups,
        source_length_unit=args.source_length_unit,
        requested=args.z_placement,
    )

    group, raw_layer_points, raw_layer_fields, raw_layer_z = select_example_layer(
        groups,
        example_z_mm=args.example_z,
        source_length_unit=args.source_length_unit,
        z_placement=resolved_z_mode,
        base_z_mm=base_z_mm,
    )

    mapper_iz = int(np.argmin(np.abs(z_axis - raw_layer_z)))
    mapper_layer_z = float(z_axis[mapper_iz])

    raw_mag = np.linalg.norm(raw_layer_fields, axis=1)
    mapper_slice_mag = np.linalg.norm(cube[:, :, mapper_iz, :], axis=2)

    # One common field scale for all four validation figures.
    vmin = 0.0
    vmax = float(max(np.max(raw_mag), np.max(mapper_slice_mag)))

    print("Creating 1/4 quarter quiver comparison ...")
    make_quarter_quiver_comparison(
        output_dir / "01_quarter_quiver_raw_vs_mapper.png",
        raw_layer_points,
        raw_layer_fields,
        cube,
        x_axis,
        y_axis,
        z_axis,
        raw_layer_z,
        quiver_step=args.quiver_step,
        raw_quiver_max=args.raw_quiver_max,
        vmin=vmin,
        vmax=vmax,
    )

    print("Creating 2/4 quarter heatmap comparison ...")
    make_quarter_heatmap_comparison(
        output_dir / "02_quarter_heatmap_raw_vs_mapper.png",
        raw_layer_points,
        raw_layer_fields,
        cube,
        x_axis,
        y_axis,
        z_axis,
        raw_layer_z,
        vmin=vmin,
        vmax=vmax,
    )

    print("Creating 3/4 full-geometry quiver comparison ...")
    make_full_quiver_comparison(
        output_dir / "03_full_quiver_raw_vs_mapper.png",
        raw_layer_points,
        raw_layer_fields,
        cube,
        x_axis,
        y_axis,
        z_axis,
        raw_layer_z,
        quiver_step=args.quiver_step,
        raw_quiver_max=max(args.raw_quiver_max, 2200),
        vmin=vmin,
        vmax=vmax,
    )

    print("Creating 4/4 full-geometry heatmap comparison ...")
    make_full_heatmap_comparison(
        output_dir / "04_full_heatmap_raw_vs_mapper.png",
        raw_layer_points,
        raw_layer_fields,
        cube,
        x_axis,
        y_axis,
        z_axis,
        raw_layer_z,
        vmin=vmin,
        vmax=vmax,
    )

    print("Calculating raw-vs-mapper numerical agreement at matching coordinates ...")
    validation = collect_numerical_validation(
        groups,
        cube,
        x_axis,
        y_axis,
        z_axis,
        source_length_unit=args.source_length_unit,
        z_placement=resolved_z_mode,
        base_z_mm=base_z_mm,
    )

    component_rows = build_component_rows(validation)
    region_rows = build_region_rows(
        validation,
        region1_end_mm=args.region1_end,
        region2_end_mm=args.region2_end,
        region3_end_mm=args.region3_end,
    )

    component_csv = output_dir / "numerical_validation_components.csv"
    region_csv = output_dir / "numerical_validation_regions.csv"
    component_png = output_dir / "05_numerical_validation_components_table.png"
    region_png = output_dir / "06_numerical_validation_regions_table.png"

    write_component_csv(component_csv, component_rows)
    write_region_csv(region_csv, region_rows)
    make_component_table_png(component_png, component_rows)
    make_region_table_png(region_png, region_rows)

    print("\nRaw .fld vs Mapper.txt validation")
    print(f"  selected physical plate : {group.label}")
    print(f"  original source z layer : {raw_layer_z:g} mm")
    print(f"  matching mapper z layer : {mapper_layer_z:g} mm")
    print(f"  source z placement      : {resolved_z_mode}")
    print(f"  mapper grid shape       : {header.shape}")
    print(
        f"  mapper spacing [mm]     : "
        f"({header.dx_mm:g}, {header.dy_mm:g}, {header.dz_mm:g})"
    )
    print(f"  source points compared  : {len(validation.points):,}")
    print("  validation B threshold  : none")

    print_component_table(component_rows)
    print_region_table(region_rows)

    print("\nSaved validation outputs")
    filenames = [
        "01_quarter_quiver_raw_vs_mapper.png",
        "02_quarter_heatmap_raw_vs_mapper.png",
        "03_full_quiver_raw_vs_mapper.png",
        "04_full_heatmap_raw_vs_mapper.png",
        "05_numerical_validation_components_table.png",
        "06_numerical_validation_regions_table.png",
        "numerical_validation_components.csv",
        "numerical_validation_regions.csv",
    ]

    for filename in filenames:
        print(f"  {output_dir / filename}")

    print("\nNotes")
    print(
        "  Raw heatmaps are binned only for display. Blank raw cells mean that "
        "no original .fld sample landed in that mapper XY cell."
    )
    print(
        "  Numerical errors are calculated directly at the original raw source "
        "coordinates; the heatmap binning is not used in the error calculation."
    )
    print(
        "  The |B| > 0.01 T active-field definition is intentionally reserved "
        "for the separate representative-B-vs-z analysis."
    )
    print(
        "  Quiver arrow lengths are normalized for direction; colour represents |B|."
    )
    print(
        "  Current masking is source-support based, not a true steel-geometry mask."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
