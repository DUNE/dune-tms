from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class MapperHeader:
    offset_x_mm: float
    offset_y_mm: float
    offset_z_mm: float
    dx_mm: float
    dy_mm: float
    dz_mm: float
    header_line_number: int
    shape: Optional[Tuple[int, int, int]] = None

    @property
    def spacing(self) -> Tuple[float, float, float]:
        return self.dx_mm, self.dy_mm, self.dz_mm


_SHAPE_RE = re.compile(r"^#\s*shape\s+(\d+)\s+(\d+)\s+(\d+)\s*$", re.I)


def read_mapper_header(path: os.PathLike | str) -> MapperHeader:
    """Read edep-sim ArbBField header and optional shape metadata."""
    path = Path(path)
    shape = None
    with path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                m = _SHAPE_RE.match(stripped)
                if m:
                    shape = tuple(map(int, m.groups()))
                continue

            vals = np.fromstring(stripped, sep=" ")
            if vals.size != 6:
                raise ValueError(
                    f"First non-comment line in {path} must contain six values "
                    "(offset_x offset_y offset_z dx dy dz)"
                )
            if np.any(vals[3:] <= 0):
                raise ValueError(f"Grid spacing must be positive in {path}")
            return MapperHeader(
                *map(float, vals),
                header_line_number=line_number,
                shape=shape,
            )
    raise ValueError(f"No edep-sim grid header found in {path}")


def iter_mapper_rows(path: os.PathLike | str) -> Iterator[np.ndarray]:
    """Yield x y z Bx By Bz Bmag rows, skipping comments and header."""
    path = Path(path)
    header = read_mapper_header(path)
    with path.open("r") as handle:
        for _ in range(header.header_line_number):
            next(handle)
        for line_number, line in enumerate(
            handle, start=header.header_line_number + 1
        ):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            vals = np.fromstring(stripped, sep=" ")
            if vals.size != 7:
                raise ValueError(
                    f"{path}:{line_number}: expected 7 data columns, got {vals.size}"
                )
            yield vals


def load_mapper_data(path: os.PathLike | str) -> Tuple[MapperHeader, np.ndarray]:
    """Load all mapper rows into memory for plotting/interactive analysis."""
    header = read_mapper_header(path)
    data = np.asarray(list(iter_mapper_rows(path)), dtype=np.float64)
    if data.size == 0:
        data = np.empty((0, 7), dtype=np.float64)
    return header, data


def validate_mapper_file(
    path: os.PathLike | str,
    *,
    atol_coord: float = 1e-6,
    atol_bmag: float = 5e-6,
) -> Dict[str, object]:
    """Stream-validate ClarkMcGrew/edep-sim ordering and uniform spacing.

    Required data ordering is z fastest, then y, then x.
    """
    path = Path(path)
    header = read_mapper_header(path)
    dx, dy, dz = header.spacing
    ox, oy, oz = header.offset_x_mm, header.offset_y_mm, header.offset_z_mm

    count = 0
    nx_seen = ny_seen = nz_seen = 0
    max_bmag_error = 0.0
    prev = None
    unique_x = set()
    unique_y = set()
    unique_z = set()

    for row in iter_mapper_rows(path):
        x, y, z, bx, by, bz, bmag = map(float, row)
        expected_mag = math.sqrt(bx * bx + by * by + bz * bz)
        max_bmag_error = max(max_bmag_error, abs(expected_mag - bmag))

        # Every coordinate must lie exactly on the declared uniform lattice.
        for value, origin, delta, label in (
            (x, ox, dx, "x"), (y, oy, dy, "y"), (z, oz, dz, "z")
        ):
            index = round((value - origin) / delta)
            reconstructed = origin + index * delta
            if abs(value - reconstructed) > atol_coord:
                raise ValueError(
                    f"Coordinate {label}={value} is off the declared lattice "
                    f"(origin={origin}, delta={delta})"
                )

        if prev is not None:
            px, py, pz = prev
            # Legal next row:
            #   same x,y and z += dz
            # or same x, y advances and z resets to oz
            # or x advances and y,z reset.
            same_xy_znext = (
                abs(x - px) <= atol_coord and
                abs(y - py) <= atol_coord and
                abs(z - (pz + dz)) <= atol_coord
            )
            same_x_ynext = (
                abs(x - px) <= atol_coord and
                y > py + atol_coord and
                abs(z - oz) <= atol_coord
            )
            xnext = (
                x > px + atol_coord and
                abs(y - oy) <= atol_coord and
                abs(z - oz) <= atol_coord
            )
            if not (same_xy_znext or same_x_ynext or xnext):
                raise ValueError(
                    "Invalid row ordering near "
                    f"x={x}, y={y}, z={z}; edep-sim requires z-fastest, "
                    "then y, then x."
                )

        unique_x.add(round(x, 9))
        unique_y.add(round(y, 9))
        unique_z.add(round(z, 9))
        prev = (x, y, z)
        count += 1

    if count == 0:
        raise ValueError(f"No data rows found in {path}")
    if max_bmag_error > atol_bmag:
        raise ValueError(
            f"Bmag column is inconsistent with Bx,By,Bz: max error "
            f"{max_bmag_error:g} T"
        )

    nx, ny, nz = len(unique_x), len(unique_y), len(unique_z)
    expected = nx * ny * nz
    if count != expected:
        raise ValueError(
            f"Grid is incomplete: rows={count}, NX*NY*NZ={expected} "
            f"({nx}*{ny}*{nz})"
        )

    if header.shape and header.shape != (nx, ny, nz):
        raise ValueError(
            f"# shape metadata {header.shape} does not match observed "
            f"{(nx, ny, nz)}"
        )

    return {
        "path": str(path),
        "rows": count,
        "shape": (nx, ny, nz),
        "offset_mm": (ox, oy, oz),
        "spacing_mm": (dx, dy, dz),
        "max_bmag_error_t": max_bmag_error,
        "valid": True,
    }
