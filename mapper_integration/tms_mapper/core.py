from __future__ import annotations

import csv
import logging
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

LOG = logging.getLogger("tms_mapper")

_FILENAME_RE = re.compile(
    r"^Plate(?P<plate>[12])_(?P<thickness>15|40|80)(?:_(?P<index>\d+))?\.fld$"
)

LengthUnit = Literal["m", "mm"]
GeometryMode = Literal["quarter", "full"]
InterpolationMode = Literal["nearest", "idw"]
ReflectionMode = Literal["axial", "polar", "none"]
ZPlacementMode = Literal["auto", "source", "zcoord"]


TMS_COORDINATE_CONVENTION = (
    "X=transverse-left-of-+Z;Y=up-opposite-gravity;Z=beam-direction"
)


def source_to_tms_axes(
    points: np.ndarray,
    fields: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return source coordinates/components in the TMS global convention.

    Marco's .fld files are already expressed in the TMS global frame:

      X : transverse horizontal axis, with +X to the left of +Z
      Y : vertical axis, positive opposite gravity
      Z : beam direction

    The corresponding magnetic-field columns are Bx, By, Bz in those same
    axes.  The mapping is therefore intentionally the identity mapping.

    Keeping this conversion as an explicit function prevents a plotting
    convention (for example, drawing Matplotlib's third axis vertically) from
    ever being mistaken for a physical coordinate permutation in Mapper.txt.
    """
    points = np.asarray(points, dtype=np.float64)
    fields = np.asarray(fields, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if fields.shape != points.shape:
        raise ValueError("fields must have the same (N, 3) shape as points")
    return points, fields


@dataclass(frozen=True)
class RawMapInfo:
    path: Path
    stem: str
    plate_part: int
    thickness_mm: float
    index: int
    z_offset_mm: float


@dataclass
class PlateGroup:
    """One physical steel plate at one longitudinal position.

    A group normally consists of the Plate1 and Plate2 quarter-geometry
    field solutions at the same z position.
    """

    z_offset_mm: float
    thickness_mm: float
    maps: Dict[int, RawMapInfo]

    @property
    def label(self) -> str:
        return f"{self.thickness_mm:g}mm@{self.z_offset_mm:g}mm"


@dataclass
class BuildSummary:
    output: Path
    geometry: str
    grid_spacing_mm: Tuple[float, float, float]
    grid_shape: Tuple[int, int, int]
    grid_bounds_mm: Tuple[float, float, float, float, float, float]
    total_points: int
    source_maps: int
    plate_groups: int
    z_placement: str
    base_z_mm: float


def _unit_scale_to_mm(unit: LengthUnit) -> float:
    if unit == "m":
        return 1000.0
    if unit == "mm":
        return 1.0
    raise ValueError(f"Unsupported source length unit: {unit}")


def read_z_coordinates(path: os.PathLike | str) -> Dict[str, float]:
    """Read Marco's z_coord.csv mapping.

    Returns a dictionary keyed by the .fld stem, e.g. ``Plate1_15_3``.
    Values are z offsets in millimetres.
    """
    result: Dict[str, float] = {}
    path = Path(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "Name" not in reader.fieldnames:
            raise ValueError(f"{path} does not contain a 'Name' column")

        z_col = None
        for candidate in ("z[mm]", "z", "Z[mm]", "Z"):
            if candidate in reader.fieldnames:
                z_col = candidate
                break
        if z_col is None:
            raise ValueError(
                f"{path} does not contain a recognized z column. "
                f"Columns: {reader.fieldnames}"
            )

        for row in reader:
            name = (row.get("Name") or "").strip()
            if not name:
                continue
            if name in result:
                raise ValueError(f"Duplicate entry for {name!r} in {path}")
            result[name] = float(row[z_col])

    if not result:
        raise ValueError(f"No z-coordinate rows found in {path}")
    return result


def discover_plate_groups(field_dir: os.PathLike | str) -> List[PlateGroup]:
    """Discover and pair Plate1/Plate2 maps using z_coord.csv."""
    field_dir = Path(field_dir)
    z_path = field_dir / "z_coord.csv"
    if not z_path.exists():
        raise FileNotFoundError(f"Missing required file: {z_path}")
    zmap = read_z_coordinates(z_path)

    grouped: Dict[Tuple[float, float], PlateGroup] = {}
    matched = 0

    for path in sorted(field_dir.glob("*.fld")):
        match = _FILENAME_RE.match(path.name)
        if not match:
            LOG.warning("Ignoring unrecognized .fld filename: %s", path.name)
            continue

        stem = path.stem
        if stem not in zmap:
            raise ValueError(f"{path.name} is missing from {z_path.name}")

        part = int(match.group("plate"))
        thickness = float(match.group("thickness"))
        index = int(match.group("index") or 0)
        z_offset = float(zmap[stem])
        key = (z_offset, thickness)

        if key not in grouped:
            grouped[key] = PlateGroup(
                z_offset_mm=z_offset,
                thickness_mm=thickness,
                maps={},
            )
        if part in grouped[key].maps:
            raise ValueError(
                f"Duplicate Plate{part} map for z={z_offset} mm, "
                f"thickness={thickness} mm"
            )
        grouped[key].maps[part] = RawMapInfo(
            path=path,
            stem=stem,
            plate_part=part,
            thickness_mm=thickness,
            index=index,
            z_offset_mm=z_offset,
        )
        matched += 1

    if matched == 0:
        raise ValueError(f"No recognized Plate1/Plate2 .fld maps found in {field_dir}")

    groups = sorted(grouped.values(), key=lambda g: (g.z_offset_mm, g.thickness_mm))

    incomplete = [g for g in groups if set(g.maps) != {1, 2}]
    if incomplete:
        details = ", ".join(f"{g.label}: parts={sorted(g.maps)}" for g in incomplete)
        raise ValueError(f"Incomplete Plate1/Plate2 pairs: {details}")

    # Cross-check that z_coord.csv did not refer to field maps that are missing.
    missing_fld = sorted(
        name for name in zmap
        if _FILENAME_RE.match(name + ".fld") and not (field_dir / (name + ".fld")).exists()
    )
    if missing_fld:
        preview = ", ".join(missing_fld[:8])
        extra = "" if len(missing_fld) <= 8 else f" ... (+{len(missing_fld)-8})"
        raise ValueError(f"z_coord.csv references missing .fld files: {preview}{extra}")

    return groups


def read_fld(
    path: os.PathLike | str,
    *,
    source_length_unit: LengthUnit = "m",
) -> Tuple[np.ndarray, np.ndarray]:
    """Read a COMSOL-style .fld file.

    The parser deliberately ignores textual header lines and accepts every
    line containing exactly six finite numeric values:

        x y z Bx By Bz

    Coordinates are returned in millimetres and B components in tesla.
    """
    path = Path(path)
    scale = _unit_scale_to_mm(source_length_unit)
    xyz: List[np.ndarray] = []
    bvec: List[np.ndarray] = []

    with path.open("r", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                values = np.fromstring(stripped, sep=" ")
            except ValueError:
                continue
            if values.size != 6:
                continue
            if not np.all(np.isfinite(values)):
                continue
            xyz.append(values[:3])
            bvec.append(values[3:6])

    if not xyz:
        raise ValueError(
            f"No six-column numeric field rows were found in {path}. "
            "Expected x y z Bx By Bz."
        )

    points = np.asarray(xyz, dtype=np.float64)
    fields = np.asarray(bvec, dtype=np.float64)
    points *= scale
    return source_to_tms_axes(points, fields)


def reflect_axial_from_quarter(
    quarter_b: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Reflect +x,+y magnetic field values into all four quadrants.

    Magnetic field is an axial (pseudo-)vector.  Starting from the +x,+y
    source quadrant, the parity is

      x reflection: (Bx, By, Bz) -> ( Bx, -By, -Bz)
      y reflection: (Bx, By, Bz) -> (-Bx,  By, -Bz)

    Therefore for sx = sign(x), sy = sign(y):

      Bx = sy * Bx_q
      By = sx * By_q
      Bz = sx * sy * Bz_q

    Coordinates on a symmetry plane use +1 for that axis.
    """
    quarter_b = np.asarray(quarter_b, dtype=np.float64)
    sx = np.where(np.asarray(x) < 0.0, -1.0, 1.0)
    sy = np.where(np.asarray(y) < 0.0, -1.0, 1.0)
    out = quarter_b.copy()
    out[:, 0] *= sy
    out[:, 1] *= sx
    out[:, 2] *= sx * sy
    return out


def _reflect_polar_from_quarter(
    quarter_b: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Diagnostic polar-vector reflection (not the physical B default)."""
    sx = np.where(np.asarray(x) < 0.0, -1.0, 1.0)
    sy = np.where(np.asarray(y) < 0.0, -1.0, 1.0)
    out = quarter_b.copy()
    out[:, 0] *= sx
    out[:, 1] *= sy
    return out


def _raw_z_min_mm(info: RawMapInfo, source_length_unit: LengthUnit) -> float:
    points, _ = read_fld(info.path, source_length_unit=source_length_unit)
    return float(np.min(points[:, 2]))


def determine_z_placement(
    groups: Sequence[PlateGroup],
    *,
    source_length_unit: LengthUnit = "m",
    requested: ZPlacementMode = "auto",
    tolerance_mm: float = 2.0,
) -> Tuple[Literal["source", "zcoord"], float]:
    """Determine whether .fld z coordinates are already globally positioned.

    ``z_coord.csv`` is always used to organize physical plates.  In ``auto``
    mode we compare the raw .fld z displacement of representative plates
    against the CSV displacement.  If they agree, raw source z is retained.
    Otherwise each .fld is shifted so its minimum z equals

        base_z + z_coord.csv offset.

    ``base_z`` is inferred from the first physical plate's raw minimum z.
    """
    if not groups:
        raise ValueError("No plate groups")

    first = groups[0]
    first_info = first.maps[1]
    base_raw = _raw_z_min_mm(first_info, source_length_unit)

    if requested == "source":
        return "source", base_raw
    if requested == "zcoord":
        return "zcoord", base_raw

    # Probe a few groups across the detector, not all 80 maps.
    probe_indices = sorted(set([0, len(groups) // 4, len(groups) // 2,
                                (3 * len(groups)) // 4, len(groups) - 1]))
    errors = []
    first_offset = groups[0].z_offset_mm
    for i in probe_indices:
        group = groups[i]
        raw = _raw_z_min_mm(group.maps[1], source_length_unit)
        raw_delta = raw - base_raw
        csv_delta = group.z_offset_mm - first_offset
        errors.append(abs(raw_delta - csv_delta))

    worst = max(errors) if errors else 0.0
    mode: Literal["source", "zcoord"] = "source" if worst <= tolerance_mm else "zcoord"
    LOG.info(
        "Auto z-placement: max source-vs-CSV displacement error %.3f mm -> %s",
        worst,
        mode,
    )
    return mode, base_raw


def _load_aligned_group(
    group: PlateGroup,
    *,
    source_length_unit: LengthUnit,
    z_placement: Literal["source", "zcoord"],
    base_z_mm: float,
) -> Tuple[np.ndarray, np.ndarray]:
    point_parts = []
    field_parts = []
    expected_start = base_z_mm + group.z_offset_mm

    for part in (1, 2):
        info = group.maps[part]
        points, fields = read_fld(
            info.path,
            source_length_unit=source_length_unit,
        )
        if z_placement == "zcoord":
            shift = expected_start - float(np.min(points[:, 2]))
            points[:, 2] += shift
        point_parts.append(points)
        field_parts.append(fields)

    points = np.concatenate(point_parts, axis=0)
    fields = np.concatenate(field_parts, axis=0)
    return points, fields


def _estimate_support_radius_xy(points: np.ndarray) -> float:
    """Estimate a conservative support radius from the quarter XY mesh."""
    xy = np.asarray(points[:, :2], dtype=np.float64)
    if len(xy) < 2:
        return math.inf

    # Deduplicate with sub-micron precision in mm coordinates.
    unique = np.unique(np.round(xy, decimals=6), axis=0)
    if len(unique) < 2:
        return math.inf

    # Limit the probe count but query against the full unique XY set.
    if len(unique) > 5000:
        idx = np.linspace(0, len(unique) - 1, 5000, dtype=int)
        probe = unique[idx]
    else:
        probe = unique

    tree = cKDTree(unique)
    k = min(5, len(unique))
    distances, _ = tree.query(probe, k=k)
    if k == 1:
        return math.inf
    distances = np.atleast_2d(distances)

    positive = distances[distances > 1e-8]
    if positive.size == 0:
        return math.inf

    typical = float(np.percentile(positive, 95))
    return 2.5 * typical


def _interpolate_tree(
    tree: cKDTree,
    source_fields: np.ndarray,
    query_points: np.ndarray,
    *,
    mode: InterpolationMode,
    k: int,
    idw_power: float,
    support_radius_mm: Optional[float],
    workers: int,
) -> np.ndarray:
    if query_points.size == 0:
        return np.empty((0, 3), dtype=np.float64)

    if mode == "nearest":
        dist, idx = tree.query(query_points, k=1, workers=workers)
        out = source_fields[idx].copy()
        nearest = np.asarray(dist)
    else:
        kk = max(1, min(int(k), len(source_fields)))
        dist, idx = tree.query(query_points, k=kk, workers=workers)
        if kk == 1:
            dist = dist[:, None]
            idx = idx[:, None]
        dist = np.asarray(dist, dtype=np.float64)
        idx = np.asarray(idx)
        nearest = dist[:, 0]

        out = np.empty((len(query_points), 3), dtype=np.float64)
        exact = dist[:, 0] <= 1e-12
        if np.any(exact):
            out[exact] = source_fields[idx[exact, 0]]

        non = ~exact
        if np.any(non):
            d = np.maximum(dist[non], 1e-12)
            weights = 1.0 / np.power(d, idw_power)
            vals = source_fields[idx[non]]
            out[non] = np.sum(vals * weights[:, :, None], axis=1) / np.sum(
                weights, axis=1
            )[:, None]

    if support_radius_mm is not None and math.isfinite(support_radius_mm):
        out[nearest > support_radius_mm] = 0.0

    return out


def _snap_min(value: float, spacing: float) -> float:
    return math.floor(value / spacing) * spacing


def _snap_max(value: float, spacing: float) -> float:
    return math.ceil(value / spacing) * spacing


def _axis_from_bounds(vmin: float, vmax: float, spacing: float) -> np.ndarray:
    if spacing <= 0:
        raise ValueError("Grid spacing must be positive")
    if vmax < vmin:
        raise ValueError(f"Invalid bounds: {vmin} > {vmax}")
    n = int(math.ceil((vmax - vmin) / spacing - 1e-12)) + 1
    return vmin + np.arange(n, dtype=np.float64) * spacing


def _resolve_bounds(
    source_bounds: Tuple[float, float, float, float, float, float],
    *,
    geometry: GeometryMode,
    dx: float,
    dy: float,
    dz: float,
    user_bounds: Optional[Tuple[Optional[float], Optional[float], Optional[float],
                                Optional[float], Optional[float], Optional[float]]] = None,
) -> Tuple[float, float, float, float, float, float]:
    qxmin, qxmax, qymin, qymax, zmin, zmax = source_bounds

    if geometry == "full":
        source_xmin, source_xmax = -max(abs(qxmin), abs(qxmax)), max(abs(qxmin), abs(qxmax))
        source_ymin, source_ymax = -max(abs(qymin), abs(qymax)), max(abs(qymin), abs(qymax))
    else:
        source_xmin, source_xmax = max(0.0, qxmin), qxmax
        source_ymin, source_ymax = max(0.0, qymin), qymax

    defaults = (
        _snap_min(source_xmin, dx),
        _snap_max(source_xmax, dx),
        _snap_min(source_ymin, dy),
        _snap_max(source_ymax, dy),
        _snap_min(zmin, dz),
        _snap_max(zmax, dz),
    )

    if user_bounds is None:
        return defaults

    resolved = []
    spacings = (dx, dx, dy, dy, dz, dz)
    is_min = (True, False, True, False, True, False)
    for user, default, spacing, minimum in zip(user_bounds, defaults, spacings, is_min):
        if user is None:
            resolved.append(default)
        else:
            # Expand outward to preserve the exact requested spacing.
            resolved.append(_snap_min(user, spacing) if minimum else _snap_max(user, spacing))
    return tuple(resolved)  # type: ignore[return-value]


def _source_bounds(
    groups: Sequence[PlateGroup],
    *,
    source_length_unit: LengthUnit,
    z_placement: Literal["source", "zcoord"],
    base_z_mm: float,
) -> Tuple[float, float, float, float, float, float]:
    xmin = ymin = zmin = math.inf
    xmax = ymax = zmax = -math.inf

    # We need the true transverse limits.  Parsing 149 MB once is cheap
    # compared with the actual uniform-grid interpolation, and avoids hard
    # coding geometry dimensions.
    for i, group in enumerate(groups, 1):
        points, _ = _load_aligned_group(
            group,
            source_length_unit=source_length_unit,
            z_placement=z_placement,
            base_z_mm=base_z_mm,
        )
        xmin = min(xmin, float(points[:, 0].min()))
        xmax = max(xmax, float(points[:, 0].max()))
        ymin = min(ymin, float(points[:, 1].min()))
        ymax = max(ymax, float(points[:, 1].max()))
        zmin = min(zmin, float(points[:, 2].min()))
        zmax = max(zmax, float(points[:, 2].max()))
        if i % 10 == 0 or i == len(groups):
            LOG.info("Scanned source bounds: %d/%d plate groups", i, len(groups))

    if not all(np.isfinite([xmin, xmax, ymin, ymax, zmin, zmax])):
        raise ValueError("Could not determine finite source bounds")
    return xmin, xmax, ymin, ymax, zmin, zmax


def _parse_support_radius(value: str | float | None, points: np.ndarray) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (float, int)):
        if float(value) <= 0:
            raise ValueError("support radius must be positive")
        return float(value)

    text = str(value).strip().lower()
    if text == "none":
        return None
    if text == "auto":
        radius = _estimate_support_radius_xy(points)
        LOG.info("Auto support radius: %.3f mm", radius)
        return radius

    radius = float(text)
    if radius <= 0:
        raise ValueError("support radius must be positive")
    return radius


def _write_edep_file(
    output: Path,
    field: np.memmap,
    x_values: np.ndarray,
    y_values: np.ndarray,
    z_values: np.ndarray,
    *,
    dx: float,
    dy: float,
    dz: float,
    metadata: Sequence[str],
    x_block: int = 2,
) -> None:
    nx, ny, nz, _ = field.shape
    with output.open("w") as handle:
        handle.write("# TMS magnetic-field map for ClarkMcGrew/edep-sim ArbBField\n")
        for line in metadata:
            handle.write(f"# {line}\n")
        handle.write(f"# shape {nx} {ny} {nz}\n")
        handle.write("# units position=mm field=tesla\n")
        handle.write("# data_columns x y z Bx By Bz Bmag\n")
        # First non-comment line: offset_x offset_y offset_z hx hy hz
        handle.write(
            f"{x_values[0]:.10g} {y_values[0]:.10g} {z_values[0]:.10g} "
            f"{dx:.10g} {dy:.10g} {dz:.10g}\n"
        )

        fmt = "%.10g %.10g %.10g %.10g %.10g %.10g %.10g"
        for ix0 in range(0, nx, x_block):
            ix1 = min(nx, ix0 + x_block)
            bx = ix1 - ix0
            b = np.asarray(field[ix0:ix1], dtype=np.float64).reshape(-1, 3)
            mag = np.linalg.norm(b, axis=1)
            xcol = np.repeat(x_values[ix0:ix1], ny * nz)
            ycol = np.tile(np.repeat(y_values, nz), bx)
            zcol = np.tile(z_values, bx * ny)
            rows = np.column_stack((xcol, ycol, zcol, b, mag))
            np.savetxt(handle, rows, fmt=fmt)


def build_uniform_mapper(
    *,
    field_dir: os.PathLike | str,
    output: os.PathLike | str,
    grid_x_size_mm: float,
    grid_y_size_mm: float,
    grid_z_size_mm: float,
    geometry: GeometryMode = "full",
    source_length_unit: LengthUnit = "m",
    z_placement: ZPlacementMode = "auto",
    z_auto_tolerance_mm: float = 2.0,
    interpolation: InterpolationMode = "idw",
    k_neighbors: int = 8,
    idw_power: float = 2.0,
    support_radius_mm: str | float | None = "auto",
    reflection: ReflectionMode = "axial",
    bounds: Optional[Tuple[Optional[float], Optional[float], Optional[float],
                           Optional[float], Optional[float], Optional[float]]] = None,
    max_points: int = 50_000_000,
    force: bool = False,
    workers: int = -1,
    chunk_size: int = 200_000,
    dry_run: bool = False,
) -> BuildSummary:
    """Build a uniformly spaced full- or quarter-TMS edep-sim field map."""
    field_dir = Path(field_dir)
    output = Path(output)
    dx, dy, dz = map(float, (grid_x_size_mm, grid_y_size_mm, grid_z_size_mm))
    if min(dx, dy, dz) <= 0:
        raise ValueError("GRID_X_SIZE, GRID_Y_SIZE and GRID_Z_SIZE must all be > 0")

    groups = discover_plate_groups(field_dir)
    resolved_z_mode, base_z_mm = determine_z_placement(
        groups,
        source_length_unit=source_length_unit,
        requested=z_placement,
        tolerance_mm=z_auto_tolerance_mm,
    )

    src_bounds = _source_bounds(
        groups,
        source_length_unit=source_length_unit,
        z_placement=resolved_z_mode,
        base_z_mm=base_z_mm,
    )
    xmin, xmax, ymin, ymax, zmin, zmax = _resolve_bounds(
        src_bounds,
        geometry=geometry,
        dx=dx,
        dy=dy,
        dz=dz,
        user_bounds=bounds,
    )

    x_values = _axis_from_bounds(xmin, xmax, dx)
    y_values = _axis_from_bounds(ymin, ymax, dy)
    z_values = _axis_from_bounds(zmin, zmax, dz)
    nx, ny, nz = len(x_values), len(y_values), len(z_values)
    total = nx * ny * nz

    summary = BuildSummary(
        output=output,
        geometry=geometry,
        grid_spacing_mm=(dx, dy, dz),
        grid_shape=(nx, ny, nz),
        grid_bounds_mm=(float(x_values[0]), float(x_values[-1]),
                        float(y_values[0]), float(y_values[-1]),
                        float(z_values[0]), float(z_values[-1])),
        total_points=total,
        source_maps=2 * len(groups),
        plate_groups=len(groups),
        z_placement=resolved_z_mode,
        base_z_mm=base_z_mm,
    )

    LOG.info("Source quarter bounds [mm]: x=(%.3f, %.3f), y=(%.3f, %.3f), z=(%.3f, %.3f)",
             *src_bounds)
    LOG.info("Output grid bounds [mm]: x=(%.3f, %.3f), y=(%.3f, %.3f), z=(%.3f, %.3f)",
             x_values[0], x_values[-1], y_values[0], y_values[-1],
             z_values[0], z_values[-1])
    LOG.info("Output grid: NX=%d NY=%d NZ=%d -> %d points", nx, ny, nz, total)

    # ~75-120 text bytes/point depending on precision and signs.
    LOG.info("Approximate text output size: %.2f GB", total * 95 / 1e9)

    if total > max_points and not force:
        raise RuntimeError(
            f"Requested grid contains {total:,} points, exceeding --max-points "
            f"{max_points:,}. Increase spacing, raise --max-points, or use --force."
        )

    if dry_run:
        return summary

    output.parent.mkdir(parents=True, exist_ok=True)
    temp = tempfile.NamedTemporaryFile(
        prefix="tms_mapper_", suffix=".bin", dir=str(output.parent), delete=False
    )
    temp_path = Path(temp.name)
    temp.close()

    field = np.memmap(
        temp_path,
        mode="w+",
        dtype=np.float32,
        shape=(nx, ny, nz, 3),
    )
    field[:] = 0.0

    try:
        xx, yy = np.meshgrid(x_values, y_values, indexing="ij")
        plane_x = xx.ravel()
        plane_y = yy.ravel()

        if geometry == "full":
            quarter_x = np.abs(plane_x)
            quarter_y = np.abs(plane_y)
        else:
            quarter_x = plane_x
            quarter_y = plane_y

        for gi, group in enumerate(groups, 1):
            points, bsrc = _load_aligned_group(
                group,
                source_length_unit=source_length_unit,
                z_placement=resolved_z_mode,
                base_z_mm=base_z_mm,
            )
            group_zmin = float(points[:, 2].min())
            group_zmax = float(points[:, 2].max())

            # Only grid planes physically inside this source plate slab.
            z_indices = np.where(
                (z_values >= group_zmin - 1e-6) &
                (z_values <= group_zmax + 1e-6)
            )[0]
            if len(z_indices) == 0:
                LOG.warning(
                    "%s has no z grid plane at dz=%.3f mm; it will not appear in output",
                    group.label, dz,
                )
                continue

            radius = _parse_support_radius(support_radius_mm, points)
            tree = cKDTree(points)

            LOG.info(
                "[%d/%d] %s: source=%d, z=[%.3f,%.3f], grid_z_planes=%d, support=%s",
                gi, len(groups), group.label, len(points), group_zmin, group_zmax,
                len(z_indices), "none" if radius is None else f"{radius:.2f} mm",
            )

            for iz in z_indices:
                qz = np.full_like(quarter_x, z_values[iz], dtype=np.float64)
                q = np.column_stack((quarter_x, quarter_y, qz))

                # Query in chunks to cap transient neighbor-array memory.
                vals = np.empty((len(q), 3), dtype=np.float64)
                for start in range(0, len(q), chunk_size):
                    stop = min(len(q), start + chunk_size)
                    vals[start:stop] = _interpolate_tree(
                        tree,
                        bsrc,
                        q[start:stop],
                        mode=interpolation,
                        k=k_neighbors,
                        idw_power=idw_power,
                        support_radius_mm=radius,
                        workers=workers,
                    )

                if geometry == "full":
                    if reflection == "axial":
                        vals = reflect_axial_from_quarter(vals, plane_x, plane_y)
                    elif reflection == "polar":
                        vals = _reflect_polar_from_quarter(vals, plane_x, plane_y)
                    elif reflection == "none":
                        pass
                    else:
                        raise ValueError(f"Unknown reflection mode: {reflection}")

                # (x,y) mesh is flattened in C order, matching reshape below.
                field[:, :, iz, :] = vals.reshape(nx, ny, 3).astype(np.float32)

            del tree, points, bsrc

        field.flush()

        metadata = [
            f"geometry={geometry}",
            f"source_dir={field_dir}",
            f"source_maps={2 * len(groups)} physical_plate_positions={len(groups)}",
            f"z_placement={resolved_z_mode} base_z_mm={base_z_mm:.10g}",
            f"interpolation={interpolation} k={k_neighbors} idw_power={idw_power}",
            f"support_radius_mm={support_radius_mm}",
            f"reflection={reflection}",
            f"coordinate_convention={TMS_COORDINATE_CONVENTION}",
            "source_axis_mapping=x->X,y->Y,z->Z;Bx->Bx,By->By,Bz->Bz",
            "ordering=z-fastest,y-next,x-slowest",
        ]
        _write_edep_file(
            output,
            field,
            x_values,
            y_values,
            z_values,
            dx=dx,
            dy=dy,
            dz=dz,
            metadata=metadata,
        )
    finally:
        try:
            del field
        except Exception:
            pass
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass

    return summary
