#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from tms_mapper.edep import load_mapper_data


_COMPONENT_INDEX = {
    "bx": 3,
    "by": 4,
    "bz": 5,
    "mag": 6,
}

_UNIT_SCALE = {
    "mm": 1.0,
    "m": 1.0e-3,
}


def _component(data: np.ndarray, name: str) -> np.ndarray:
    return data[:, _COMPONENT_INDEX[name]]


def _component_label(name: str) -> str:
    return "|B| [T]" if name == "mag" else f"B{name[-1]} [T]"


def _unit_scale(units: str) -> float:
    return _UNIT_SCALE[units]


def _slice_request_to_mm(value: float | None, units: str) -> float | None:
    if value is None:
        return None
    return float(value) / _unit_scale(units)


def _quarter_filter(data: np.ndarray) -> np.ndarray:
    eps = 1e-9
    return data[(data[:, 0] >= -eps) & (data[:, 1] >= -eps)]


def _expand_quarter_for_plot(data: np.ndarray) -> np.ndarray:
    """Reflect quarter data into four quadrants with axial-vector parity."""
    pieces = []
    for sx, sy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        d = data.copy()
        d[:, 0] *= sx
        d[:, 1] *= sy
        d[:, 3] *= sy       # Bx
        d[:, 4] *= sx       # By
        d[:, 5] *= sx * sy  # Bz
        d[:, 6] = np.linalg.norm(d[:, 3:6], axis=1)
        pieces.append(d)

    out = np.concatenate(pieces, axis=0)
    key = np.round(out[:, :3], 9)
    _, idx = np.unique(key, axis=0, return_index=True)
    return out[np.sort(idx)]


def _apply_geometry(data: np.ndarray, geometry: str) -> np.ndarray:
    if geometry == "stored":
        return data
    if geometry == "quarter":
        return _quarter_filter(data)
    if geometry == "full":
        if np.min(data[:, 0]) >= -1e-9 and np.min(data[:, 1]) >= -1e-9:
            return _expand_quarter_for_plot(data)
        return data
    raise ValueError(geometry)


def _apply_ranges(data: np.ndarray, args) -> np.ndarray:
    """Crop to optional X/Y/Z ranges supplied in the requested display units."""
    scale = _unit_scale(args.units)
    mask = np.ones(len(data), dtype=bool)
    for axis_col, values in enumerate((args.x_range, args.y_range, args.z_range)):
        if values is None:
            continue
        lo_display, hi_display = map(float, values)
        if hi_display < lo_display:
            lo_display, hi_display = hi_display, lo_display
        lo_mm = lo_display / scale
        hi_mm = hi_display / scale
        mask &= (data[:, axis_col] >= lo_mm) & (data[:, axis_col] <= hi_mm)
    return data[mask]


def _nearest_slice_value(data: np.ndarray, axis_col: int, requested: float | None) -> float:
    vals = np.unique(data[:, axis_col])
    if requested is None:
        requested = float(vals[len(vals) // 2])
    return float(vals[np.argmin(np.abs(vals - requested))])


def _slice_for_plane(
    data: np.ndarray,
    plane: str,
    requested_mm: float | None,
) -> Tuple[np.ndarray, float, Tuple[int, int], str]:
    """Return one regular-grid slice in the physical TMS display convention.

    +Z is the beam direction, +Y is up/opposite gravity, and +X is transverse.
    For longitudinal 2D views Z is drawn horizontally.
    """
    if plane == "xy":
        fixed_col, dims, fixed_axis = 2, (0, 1), "Z"
    elif plane == "xz":
        fixed_col, dims, fixed_axis = 1, (2, 0), "Y"
    elif plane == "yz":
        fixed_col, dims, fixed_axis = 0, (2, 1), "X"
    else:
        raise ValueError(plane)

    fixed = _nearest_slice_value(data, fixed_col, requested_mm)
    mask = np.isclose(data[:, fixed_col], fixed, atol=1e-8)
    return data[mask], fixed, dims, fixed_axis


def _axis_label(position_col: int, units: str) -> str:
    labels = {
        0: f"X [{units}] — transverse",
        1: f"Y [{units}] — up / -gravity",
        2: f"Z [{units}] — beam",
    }
    return labels[position_col]


def _grid_matrix(slice_data: np.ndarray, a_col: int, b_col: int, values: np.ndarray):
    a = np.unique(slice_data[:, a_col])
    b = np.unique(slice_data[:, b_col])
    matrix = np.full((len(b), len(a)), np.nan, dtype=float)
    amap = {round(v, 9): i for i, v in enumerate(a)}
    bmap = {round(v, 9): i for i, v in enumerate(b)}
    for row, value in zip(slice_data, values):
        matrix[bmap[round(row[b_col], 9)], amap[round(row[a_col], 9)]] = value
    return a, b, matrix


def _default_title(args, kind: str, fixed_axis: str | None = None, fixed_display: float | None = None) -> str:
    if args.title:
        return args.title
    component = "|B|" if args.component == "mag" else f"B{args.component[-1]}"
    if kind == "3d":
        return (
            f"TMS 3D magnetic field sample — colored by {component}\n"
            "Z = beam, Y = up (-gravity), X = transverse"
        )
    if kind == "quiver":
        return f"Magnetic Field Vector Plot ({args.plane.upper()} Plane at {fixed_axis} ≈ {fixed_display:g} {args.units})"
    if kind == "scatter":
        return f"{component} Across Uniform Grid ({args.plane.upper()} Plane at {fixed_axis} ≈ {fixed_display:g} {args.units})"
    return f"TMS {component} — {args.plane.upper()} slice at {fixed_axis} ≈ {fixed_display:g} {args.units}"


def _finish_2d_axes(ax, args, dims: Tuple[int, int]) -> None:
    ax.set_xlabel(_axis_label(dims[0], args.units))
    ax.set_ylabel(_axis_label(dims[1], args.units))
    ax.set_aspect("equal")
    if not args.no_grid:
        ax.grid(True, alpha=0.28, linestyle="--", linewidth=0.6)


def plot_2d(args, data: np.ndarray):
    requested_mm = _slice_request_to_mm(args.slice, args.units)
    s, fixed_mm, dims, fixed_axis = _slice_for_plane(data, args.plane, requested_mm)
    values = _component(s, args.component)
    a, b, matrix = _grid_matrix(s, dims[0], dims[1], values)
    scale = _unit_scale(args.units)
    a = a * scale
    b = b * scale
    fixed_display = fixed_mm * scale

    fig, ax = plt.subplots(figsize=(9, 7))
    mesh = ax.pcolormesh(
        a, b, matrix,
        shading="auto",
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
    )
    fig.colorbar(mesh, ax=ax, label=_component_label(args.component))
    _finish_2d_axes(ax, args, dims)
    ax.set_title(_default_title(args, "2d", fixed_axis, fixed_display))
    fig.tight_layout()
    return fig


def plot_scatter(args, data: np.ndarray):
    """Scatter regular-grid coordinates, colored by field strength/component."""
    requested_mm = _slice_request_to_mm(args.slice, args.units)
    s, fixed_mm, dims, fixed_axis = _slice_for_plane(data, args.plane, requested_mm)
    scale = _unit_scale(args.units)
    fixed_display = fixed_mm * scale
    values = _component(s, args.component)

    fig, ax = plt.subplots(figsize=(9, 7))
    sc = ax.scatter(
        s[:, dims[0]] * scale,
        s[:, dims[1]] * scale,
        c=values,
        s=args.scatter_size,
        alpha=args.alpha,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        linewidths=0,
    )
    fig.colorbar(sc, ax=ax, label=_component_label(args.component))
    _finish_2d_axes(ax, args, dims)
    ax.set_title(_default_title(args, "scatter", fixed_axis, fixed_display))
    fig.tight_layout()
    return fig


def plot_quiver(args, data: np.ndarray):
    requested_mm = _slice_request_to_mm(args.slice, args.units)
    s, fixed_mm, dims, fixed_axis = _slice_for_plane(data, args.plane, requested_mm)
    stride = max(1, args.stride)

    a_vals = np.unique(s[:, dims[0]])[::stride]
    b_vals = np.unique(s[:, dims[1]])[::stride]
    keep_a = np.isin(s[:, dims[0]], a_vals)
    keep_b = np.isin(s[:, dims[1]], b_vals)
    s = s[keep_a & keep_b]

    field_col_for_position_col = {0: 3, 1: 4, 2: 5}
    u = s[:, field_col_for_position_col[dims[0]]].astype(float, copy=True)
    v = s[:, field_col_for_position_col[dims[1]]].astype(float, copy=True)

    if args.normalize_arrows:
        norm = np.hypot(u, v)
        nonzero = norm > 0
        u[nonzero] /= norm[nonzero]
        v[nonzero] /= norm[nonzero]

    colors = _component(s, args.component)
    scale = _unit_scale(args.units)
    fixed_display = fixed_mm * scale

    fig, ax = plt.subplots(figsize=(9, 7))
    q = ax.quiver(
        s[:, dims[0]] * scale,
        s[:, dims[1]] * scale,
        u,
        v,
        colors,
        cmap=args.cmap,
        angles="xy",
        scale_units="xy",
        scale=args.quiver_scale,
        pivot="mid",
        width=args.quiver_width,
    )
    if args.vmin is not None or args.vmax is not None:
        q.set_clim(args.vmin, args.vmax)
    fig.colorbar(q, ax=ax, label=_component_label(args.component))
    _finish_2d_axes(ax, args, dims)
    ax.set_title(_default_title(args, "quiver", fixed_axis, fixed_display))
    fig.tight_layout()
    return fig


def _xyz_for_3d(data: np.ndarray, units: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return plotting axes as X horizontal, Z beam/depth, Y vertical."""
    scale = _unit_scale(units)
    return data[:, 0] * scale, data[:, 2] * scale, data[:, 1] * scale


def plot_3d(args, data: np.ndarray):
    stride = max(1, args.stride)
    d = data[::stride]
    if args.min_field > 0:
        d = d[d[:, 6] >= args.min_field]
    if len(d) > args.max_plot_points:
        idx = np.linspace(0, len(d) - 1, args.max_plot_points, dtype=int)
        d = d[idx]
    if len(d) == 0:
        raise ValueError("No points remain after 3D filters")

    values = _component(d, args.component)
    x_plot, z_plot, y_plot = _xyz_for_3d(d, args.units)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(
        x_plot,
        z_plot,
        y_plot,
        c=values,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        s=args.point_size,
        alpha=args.alpha,
        linewidths=0,
    )
    fig.colorbar(sc, ax=ax, shrink=0.7, label=_component_label(args.component))
    ax.set_xlabel(f"X [{args.units}] — transverse")
    ax.set_ylabel(f"Z [{args.units}] — beam")
    ax.set_zlabel(f"Y [{args.units}] — up / -gravity")

    x_span = max(float(np.ptp(x_plot)), 1e-12)
    z_span = max(float(np.ptp(z_plot)), 1e-12)
    y_span = max(float(np.ptp(y_plot)), 1e-12)
    ax.set_box_aspect((x_span, z_span, y_span))
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_title(_default_title(args, "3d"))
    fig.tight_layout()
    return fig


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "CLI plotting utility for TMS Mapper.txt files. Supports 2D heatmaps, "
            "uniform-grid scatter plots, vector quiver plots, and 3D scatter plots. "
            "Plots use TMS axes: +Z beam, +Y opposite gravity, +X transverse."
        )
    )
    p.add_argument("mapper", help="Mapper.txt generated by build_mapper.py")
    p.add_argument("--kind", choices=("2d", "scatter", "quiver", "3d"), required=True)
    p.add_argument("--geometry", choices=("stored", "quarter", "full"), default="stored")
    p.add_argument("--component", choices=("mag", "bx", "by", "bz"), default="mag")
    p.add_argument("--plane", choices=("xy", "xz", "yz"), default="xy",
                   help="Slice plane for 2d/scatter/quiver")
    p.add_argument("--slice", type=float, default=None,
                   help="Orthogonal slice coordinate in --units; nearest grid plane is used")
    p.add_argument("--units", choices=("mm", "m"), default="mm",
                   help="Display units for coordinates and --slice (default: mm)")

    p.add_argument("--x-range", nargs=2, type=float, metavar=("MIN", "MAX"),
                   help="Optional X crop in --units")
    p.add_argument("--y-range", nargs=2, type=float, metavar=("MIN", "MAX"),
                   help="Optional Y crop in --units")
    p.add_argument("--z-range", nargs=2, type=float, metavar=("MIN", "MAX"),
                   help="Optional Z crop in --units")

    p.add_argument("--stride", type=int, default=3,
                   help="Grid thinning for quiver/3D (default: 3)")
    p.add_argument("--scatter-size", type=float, default=14.0,
                   help="Marker size for --kind scatter")
    p.add_argument("--point-size", type=float, default=2.0,
                   help="Marker size for --kind 3d")
    p.add_argument("--alpha", type=float, default=0.75)
    p.add_argument("--cmap", default="viridis")
    p.add_argument("--vmin", type=float, default=None,
                   help="Optional color scale minimum in tesla")
    p.add_argument("--vmax", type=float, default=None,
                   help="Optional color scale maximum in tesla")

    p.add_argument("--normalize-arrows", action="store_true",
                   help="Normalize quiver arrow lengths; color still carries field strength")
    p.add_argument("--quiver-scale", type=float, default=None,
                   help="Matplotlib quiver scale; default lets matplotlib choose")
    p.add_argument("--quiver-width", type=float, default=0.0025)

    p.add_argument("--min-field", type=float, default=0.0,
                   help="3D: omit points with |B| below this value in tesla")
    p.add_argument("--max-plot-points", type=int, default=200_000,
                   help="3D rendering cap after stride/filter")
    p.add_argument("--elev", type=float, default=22.0,
                   help="3D camera elevation in degrees")
    p.add_argument("--azim", type=float, default=-58.0,
                   help="3D camera azimuth in degrees")

    p.add_argument("--title", default=None, help="Optional custom plot title")
    p.add_argument("--no-grid", action="store_true", help="Disable 2D/quiver/scatter grid lines")
    p.add_argument("-o", "--output", default=None,
                   help="PNG/PDF/SVG output. If omitted, display interactively.")
    p.add_argument("--dpi", type=int, default=180)
    return p


def main() -> int:
    args = make_parser().parse_args()
    _, data = load_mapper_data(args.mapper)
    if len(data) == 0:
        raise SystemExit("Mapper contains no data rows")

    data = _apply_geometry(data, args.geometry)
    data = _apply_ranges(data, args)
    if len(data) == 0:
        raise SystemExit("No mapper points remain after geometry/range filters")

    if args.kind == "2d":
        fig = plot_2d(args, data)
    elif args.kind == "scatter":
        fig = plot_scatter(args, data)
    elif args.kind == "quiver":
        fig = plot_quiver(args, data)
    else:
        fig = plot_3d(args, data)

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        print(f"Wrote {out}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
