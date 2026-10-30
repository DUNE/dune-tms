
from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from tms_mapper.core import discover_plate_groups, read_fld


DUNENDGGD_REPOSITORY = "https://github.com/DUNE/dunendggd.git"
DUNENDGGD_COMMIT = "071b712939697093d4c811c5d80dbf6359b0d9d3"

DEFAULT_FIELD_VOLUMES = (
    "thinvolTMS",
    "thinvol2TMS",
    "thickvolTMS",
    "thickvol2TMS",
    "doublevolTMS",
    "doublevol2TMS",
)

_TARGET_OUTPUTS = {
    "tms_nosand": "nd_hall_with_lar_tms_nosand.gdml",
    "tms": "nd_hall_with_lar_tms_sand_stt1.gdml",
    "tms_drift1": "nd_hall_with_lar_tms_sand_drift1.gdml",
}


@dataclass(frozen=True)
class MapperHeader:
    offset_mm: Tuple[float, float, float]
    spacing_mm: Tuple[float, float, float]
    line_number: int


@dataclass
class TMSGeometryInfo:
    gdml_path: Path
    tms_volume: str
    global_translation_mm: Tuple[float, float, float]
    global_rotation: np.ndarray
    layer_centers_mm: Dict[float, np.ndarray]

    @property
    def axis_aligned(self) -> bool:
        return bool(np.allclose(self.global_rotation, np.eye(3), atol=1e-9))


@dataclass
class AlignmentReport:
    field_to_tms_local_mm: Tuple[float, float, float]
    field_to_global_mm: Tuple[float, float, float]
    max_z_residual_mm: float
    rms_z_residual_mm: float
    family_counts: Dict[float, Tuple[int, int]]
    field_pitch_mm: Dict[float, Optional[float]]
    gdml_pitch_mm: Dict[float, Optional[float]]
    geometry_compatible: bool
    messages: List[str]


def _local_name(tag: str) -> str:
    return tag.split("}", 1)[-1]


def _children(element: ET.Element, name: str) -> Iterable[ET.Element]:
    for child in element:
        if _local_name(child.tag) == name:
            yield child


def _first_child(element: ET.Element, name: str) -> Optional[ET.Element]:
    return next(iter(_children(element, name)), None)


def _length_scale_to_mm(unit: Optional[str]) -> float:
    unit = (unit or "mm").strip().lower()
    scales = {
        "mm": 1.0,
        "millimeter": 1.0,
        "millimetre": 1.0,
        "cm": 10.0,
        "centimeter": 10.0,
        "centimetre": 10.0,
        "m": 1000.0,
        "meter": 1000.0,
        "metre": 1000.0,
    }
    if unit not in scales:
        raise ValueError(f"Unsupported GDML length unit {unit!r}")
    return scales[unit]


def _angle_scale_to_rad(unit: Optional[str]) -> float:
    unit = (unit or "deg").strip().lower()
    if unit in {"deg", "degree", "degrees"}:
        return math.pi / 180.0
    if unit in {"rad", "radian", "radians"}:
        return 1.0
    raise ValueError(f"Unsupported GDML angle unit {unit!r}")


def _vector_from_element(element: ET.Element, *, angle: bool = False) -> np.ndarray:
    unit = element.get("unit") or element.get("lunit") or element.get("aunit")
    scale = _angle_scale_to_rad(unit) if angle else _length_scale_to_mm(unit)
    return np.asarray(
        [
            float(element.get("x", "0")) * scale,
            float(element.get("y", "0")) * scale,
            float(element.get("z", "0")) * scale,
        ],
        dtype=float,
    )


def _rotation_matrix_xyz(angles: np.ndarray) -> np.ndarray:
    x, y, z = map(float, angles)
    cx, sx = math.cos(x), math.sin(x)
    cy, sy = math.cos(y), math.sin(y)
    cz, sz = math.cos(z), math.sin(z)

    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], dtype=float)
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], dtype=float)
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], dtype=float)
    return rz @ ry @ rx


def _definitions(root: ET.Element) -> Tuple[Dict[str, ET.Element], Dict[str, ET.Element]]:
    positions: Dict[str, ET.Element] = {}
    rotations: Dict[str, ET.Element] = {}
    for element in root.iter():
        name = element.get("name")
        if not name:
            continue
        tag = _local_name(element.tag)
        if tag == "position":
            positions[name] = element
        elif tag == "rotation":
            rotations[name] = element
    return positions, rotations


def _physvol_transform(
    physvol: ET.Element,
    positions: Mapping[str, ET.Element],
    rotations: Mapping[str, ET.Element],
) -> Tuple[np.ndarray, np.ndarray]:
    position = _first_child(physvol, "position")
    if position is None:
        pref = _first_child(physvol, "positionref")
        if pref is not None:
            ref = pref.get("ref")
            if ref not in positions:
                raise ValueError(f"Undefined GDML positionref {ref!r}")
            position = positions[ref]
    p = np.zeros(3) if position is None else _vector_from_element(position)

    rotation = _first_child(physvol, "rotation")
    if rotation is None:
        rref = _first_child(physvol, "rotationref")
        if rref is not None:
            ref = rref.get("ref")
            if ref not in rotations:
                raise ValueError(f"Undefined GDML rotationref {ref!r}")
            rotation = rotations[ref]
    r = np.eye(3) if rotation is None else _rotation_matrix_xyz(
        _vector_from_element(rotation, angle=True)
    )
    return p, r


def _volume_ref(physvol: ET.Element) -> Optional[str]:
    vref = _first_child(physvol, "volumeref")
    return None if vref is None else vref.get("ref")


def _volume_table(root: ET.Element) -> Dict[str, ET.Element]:
    table: Dict[str, ET.Element] = {}
    for element in root.iter():
        if _local_name(element.tag) != "volume":
            continue
        name = element.get("name")
        if name:
            table[name] = element
    return table


def _world_volume_name(root: ET.Element) -> str:
    for setup in root.iter():
        if _local_name(setup.tag) != "setup":
            continue
        world = _first_child(setup, "world")
        if world is not None and world.get("ref"):
            return str(world.get("ref"))
    raise ValueError("GDML does not contain <setup><world ref=...>")


def _find_global_transform(
    root: ET.Element,
    target_volume: str,
) -> Tuple[np.ndarray, np.ndarray]:
    positions, rotations = _definitions(root)
    volumes = _volume_table(root)
    world = _world_volume_name(root)
    if world not in volumes:
        raise ValueError(f"World logical volume {world!r} not found")

    matches: List[Tuple[np.ndarray, np.ndarray]] = []

    def walk(
        volume_name: str,
        parent_t: np.ndarray,
        parent_r: np.ndarray,
        stack: Tuple[str, ...],
    ):
        if volume_name in stack:
            return
        volume = volumes.get(volume_name)
        if volume is None:
            return
        for physvol in _children(volume, "physvol"):
            child_name = _volume_ref(physvol)
            if not child_name:
                continue
            p_local, r_local = _physvol_transform(physvol, positions, rotations)
            child_t = parent_t + parent_r @ p_local
            child_r = parent_r @ r_local
            if child_name == target_volume:
                matches.append((child_t, child_r))
            walk(child_name, child_t, child_r, stack + (volume_name,))

    if world == target_volume:
        return np.zeros(3), np.eye(3)

    walk(world, np.zeros(3), np.eye(3), tuple())

    if not matches:
        raise ValueError(f"Could not locate logical volume {target_volume!r} from GDML world")
    if len(matches) != 1:
        raise ValueError(
            f"Expected one placement of {target_volume!r}; found {len(matches)}"
        )
    return matches[0]


def _layer_family_from_ref(ref: str) -> Optional[float]:
    lower = ref.lower()
    if "thinlayervol" in lower:
        return 15.0
    if "thicklayervol" in lower:
        return 40.0
    if "doublelayervol" in lower:
        return 80.0
    return None


def _extract_layer_centers(
    root: ET.Element,
    tms_volume: str,
) -> Dict[float, np.ndarray]:
    positions, rotations = _definitions(root)
    volumes = _volume_table(root)
    volume = volumes.get(tms_volume)
    if volume is None:
        raise ValueError(f"Logical volume {tms_volume!r} is missing")

    centers: Dict[float, List[np.ndarray]] = {15.0: [], 40.0: [], 80.0: []}
    for physvol in _children(volume, "physvol"):
        ref = _volume_ref(physvol)
        if not ref:
            continue
        family = _layer_family_from_ref(ref)
        if family is None:
            continue
        p, r = _physvol_transform(physvol, positions, rotations)
        if not np.allclose(r, np.eye(3), atol=1e-9):
            raise ValueError(
                f"Steel layer {ref!r} has a non-identity local rotation; "
                "automatic field alignment currently supports axis-aligned layers only."
            )
        centers[family].append(p)

    result: Dict[float, np.ndarray] = {}
    for thickness, values in centers.items():
        if values:
            arr = np.asarray(values, dtype=float)
            result[thickness] = arr[np.argsort(arr[:, 2])]
    return result


def inspect_tms_geometry(
    gdml_path: os.PathLike | str,
    *,
    tms_volume: str = "volTMS",
) -> TMSGeometryInfo:
    gdml_path = Path(gdml_path)
    root = ET.parse(gdml_path).getroot()
    t, r = _find_global_transform(root, tms_volume)
    layers = _extract_layer_centers(root, tms_volume)
    return TMSGeometryInfo(
        gdml_path=gdml_path,
        tms_volume=tms_volume,
        global_translation_mm=tuple(map(float, t)),
        global_rotation=r,
        layer_centers_mm=layers,
    )


def read_mapper_header(path: os.PathLike | str) -> MapperHeader:
    path = Path(path)
    with path.open("r") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = np.fromstring(stripped, sep=" ")
            if values.size != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected mapper header "
                    "offset_x offset_y offset_z dx dy dz"
                )
            return MapperHeader(
                offset_mm=tuple(map(float, values[:3])),
                spacing_mm=tuple(map(float, values[3:])),
                line_number=line_number,
            )
    raise ValueError(f"No mapper header found in {path}")


def _mapper_metadata(path: os.PathLike | str) -> Dict[str, str]:
    metadata: Dict[str, str] = {}
    token_re = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=([^ ]+)")
    with Path(path).open("r") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            for key, value in token_re.findall(line):
                metadata[key] = value
    return metadata


def _field_plate_centers(
    field_dir: os.PathLike | str,
    mapper_path: os.PathLike | str,
    *,
    source_length_unit: str = "m",
) -> Dict[float, np.ndarray]:
    groups = discover_plate_groups(field_dir)
    metadata = _mapper_metadata(mapper_path)
    mode = metadata.get("z_placement", "zcoord")
    base_text = metadata.get("base_z_mm")
    centers: Dict[float, List[np.ndarray]] = {15.0: [], 40.0: [], 80.0: []}

    if mode == "zcoord":
        if base_text is None:
            raise ValueError(
                "Mapper metadata does not contain base_z_mm. Regenerate it with "
                "the current build_mapper.py so GDML alignment is reproducible."
            )
        base_z = float(base_text)
        for group in groups:
            z = base_z + group.z_offset_mm + 0.5 * group.thickness_mm
            centers[group.thickness_mm].append(np.array([0.0, 0.0, z]))
    elif mode == "source":
        for group in groups:
            info = group.maps[1]
            points, _ = read_fld(info.path, source_length_unit=source_length_unit)
            z = 0.5 * (float(points[:, 2].min()) + float(points[:, 2].max()))
            centers[group.thickness_mm].append(np.array([0.0, 0.0, z]))
    else:
        raise ValueError(f"Unsupported mapper z_placement metadata {mode!r}")

    result: Dict[float, np.ndarray] = {}
    for thickness, values in centers.items():
        if values:
            arr = np.asarray(values, dtype=float)
            result[thickness] = arr[np.argsort(arr[:, 2])]
    return result


def _median_pitch(centers: np.ndarray) -> Optional[float]:
    if len(centers) < 2:
        return None
    return float(np.median(np.diff(np.sort(centers[:, 2]))))


def compute_field_alignment(
    gdml_path: os.PathLike | str,
    mapper_path: os.PathLike | str,
    field_dir: os.PathLike | str,
    *,
    tms_volume: str = "volTMS",
    source_length_unit: str = "m",
    tolerance_mm: float = 5.0,
) -> Tuple[TMSGeometryInfo, AlignmentReport]:
    geom = inspect_tms_geometry(gdml_path, tms_volume=tms_volume)
    if not geom.axis_aligned:
        raise ValueError(
            "The TMS is rotated in global GDML coordinates. edep-sim ArbBField "
            "grid coordinates are global and this integration currently applies "
            "translation only. Refusing to inject a silently rotated field."
        )

    field_centers = _field_plate_centers(
        field_dir,
        mapper_path,
        source_length_unit=source_length_unit,
    )

    messages: List[str] = []
    family_counts: Dict[float, Tuple[int, int]] = {}
    field_pitch: Dict[float, Optional[float]] = {}
    gdml_pitch: Dict[float, Optional[float]] = {}
    deltas: List[float] = []
    count_ok = True
    all_gdml_xy: List[np.ndarray] = []

    for thickness in (15.0, 40.0, 80.0):
        f = field_centers.get(thickness, np.empty((0, 3)))
        g = geom.layer_centers_mm.get(thickness, np.empty((0, 3)))
        family_counts[thickness] = (len(f), len(g))
        field_pitch[thickness] = _median_pitch(f)
        gdml_pitch[thickness] = _median_pitch(g)

        if len(g):
            all_gdml_xy.append(g[:, :2])

        if len(f) != len(g):
            count_ok = False
            messages.append(
                f"{thickness:g} mm family count mismatch: field={len(f)}, gdml={len(g)}"
            )
            continue

        if len(f):
            deltas.extend((g[:, 2] - f[:, 2]).tolist())

        fp, gp = field_pitch[thickness], gdml_pitch[thickness]
        if fp is not None and gp is not None and abs(fp - gp) > tolerance_mm:
            messages.append(
                f"{thickness:g} mm family pitch differs: "
                f"field={fp:.3f} mm, gdml={gp:.3f} mm"
            )

    if not deltas:
        raise ValueError("Could not match any field plate centers to GDML steel layers")

    z_shift = float(np.median(deltas))
    residuals = np.asarray(deltas, dtype=float) - z_shift
    max_residual = float(np.max(np.abs(residuals)))
    rms_residual = float(np.sqrt(np.mean(residuals**2)))

    if all_gdml_xy:
        xy = np.concatenate(all_gdml_xy, axis=0)
        x_shift = float(np.median(xy[:, 0]))
        y_shift = float(np.median(xy[:, 1]))
    else:
        x_shift = 0.0
        y_shift = 0.0

    local_shift = np.asarray([x_shift, y_shift, z_shift], dtype=float)
    global_shift = np.asarray(geom.global_translation_mm, dtype=float) + local_shift

    if max_residual > tolerance_mm:
        messages.append(
            f"single-translation Z alignment residual is {max_residual:.3f} mm "
            f"(tolerance {tolerance_mm:.3f} mm)"
        )

    compatible = count_ok and max_residual <= tolerance_mm
    for thickness in (15.0, 40.0, 80.0):
        fp, gp = field_pitch[thickness], gdml_pitch[thickness]
        if fp is not None and gp is not None and abs(fp - gp) > tolerance_mm:
            compatible = False

    if compatible:
        messages.append("field plate geometry is translation-compatible with GDML")
    else:
        messages.append(
            "field/GDML geometry is not exactly translation-compatible; do not "
            "treat the injected geometry as production-valid without review"
        )

    return geom, AlignmentReport(
        field_to_tms_local_mm=tuple(map(float, local_shift)),
        field_to_global_mm=tuple(map(float, global_shift)),
        max_z_residual_mm=max_residual,
        rms_z_residual_mm=rms_residual,
        family_counts=family_counts,
        field_pitch_mm=field_pitch,
        gdml_pitch_mm=gdml_pitch,
        geometry_compatible=compatible,
        messages=messages,
    )


def globalize_mapper(
    input_mapper: os.PathLike | str,
    output_mapper: os.PathLike | str,
    translation_mm: Sequence[float],
    *,
    overwrite: bool = False,
) -> MapperHeader:
    input_mapper = Path(input_mapper)
    output_mapper = Path(output_mapper)
    if input_mapper.resolve() == output_mapper.resolve():
        raise ValueError("Refusing to overwrite the source mapper in place")
    if output_mapper.exists() and not overwrite:
        raise FileExistsError(f"{output_mapper} already exists (use overwrite=True)")

    shift = np.asarray(translation_mm, dtype=float)
    if shift.shape != (3,):
        raise ValueError("translation_mm must contain exactly X,Y,Z")

    output_mapper.parent.mkdir(parents=True, exist_ok=True)
    header_written = False

    with input_mapper.open("r") as src, output_mapper.open("w") as dst:
        for line in src:
            stripped = line.strip()
            if not header_written and (not stripped or stripped.startswith("#")):
                dst.write(line)
                continue

            if not header_written:
                values = np.fromstring(stripped, sep=" ")
                if values.size != 6:
                    raise ValueError("Invalid mapper header while globalizing")
                values[:3] += shift
                dst.write(" ".join(f"{v:.10g}" for v in values) + "\n")
                dst.write(
                    "# gdml_global_translation_mm="
                    + ",".join(f"{v:.10g}" for v in shift)
                    + "\n"
                )
                dst.write(
                    "# mapper_frame=GDML-global; X/Y/Z coordinates and header "
                    "offsets include GDML alignment\n"
                )
                header_written = True
                continue

            if not stripped or stripped.startswith("#"):
                dst.write(line)
                continue

            values = np.fromstring(stripped, sep=" ")
            if values.size != 7:
                raise ValueError(
                    "Mapper data row must contain X Y Z Bx By Bz Bmag"
                )
            values[:3] += shift
            dst.write(" ".join(f"{v:.10g}" for v in values) + "\n")

    return read_mapper_header(output_mapper)


def _target_volume_elements(
    root: ET.Element,
    volume_names: Sequence[str],
) -> Dict[str, ET.Element]:
    wanted = set(volume_names)
    found: Dict[str, ET.Element] = {}
    for element in root.iter():
        if _local_name(element.tag) != "volume":
            continue
        name = element.get("name")
        if name in wanted:
            found[name] = element
    missing = sorted(wanted - set(found))
    if missing:
        raise ValueError(
            "GDML is missing expected TMS steel logical volumes: " + ", ".join(missing)
        )
    return found


def inject_arb_bfield(
    input_gdml: os.PathLike | str,
    output_gdml: os.PathLike | str,
    mapper_runtime_path: str,
    *,
    volume_names: Sequence[str] = DEFAULT_FIELD_VOLUMES,
    overwrite: bool = False,
) -> int:
    input_gdml = Path(input_gdml)
    output_gdml = Path(output_gdml)
    if input_gdml.resolve() == output_gdml.resolve():
        raise ValueError("Refusing to overwrite source GDML in place")
    if output_gdml.exists() and not overwrite:
        raise FileExistsError(f"{output_gdml} already exists (use overwrite=True)")

    tree = ET.parse(input_gdml)
    root = tree.getroot()
    volumes = _target_volume_elements(root, volume_names)

    for volume in volumes.values():
        for child in list(volume):
            if _local_name(child.tag) != "auxiliary":
                continue
            auxtype = (child.get("auxtype") or child.get("auxType") or "").strip()
            if auxtype in {"BField", "ArbBField"}:
                volume.remove(child)
        namespace = ""
        if volume.tag.startswith("{"):
            namespace = volume.tag.split("}", 1)[0] + "}"
        ET.SubElement(
            volume,
            namespace + "auxiliary",
            {"auxtype": "ArbBField", "auxvalue": str(mapper_runtime_path)},
        )

    output_gdml.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_gdml, encoding="utf-8", xml_declaration=True)
    return len(volumes)


def validate_injected_gdml(
    gdml_path: os.PathLike | str,
    *,
    volume_names: Sequence[str] = DEFAULT_FIELD_VOLUMES,
    mapper_runtime_path: Optional[str] = None,
) -> Dict[str, object]:
    root = ET.parse(gdml_path).getroot()
    volumes = _target_volume_elements(root, volume_names)
    details: Dict[str, str] = {}

    for name, volume in volumes.items():
        constant = []
        arbitrary = []
        for child in _children(volume, "auxiliary"):
            auxtype = (child.get("auxtype") or child.get("auxType") or "").strip()
            value = child.get("auxvalue") or child.get("auxValue") or ""
            if auxtype == "BField":
                constant.append(value)
            elif auxtype == "ArbBField":
                arbitrary.append(value)

        if constant:
            raise ValueError(f"{name} still contains constant BField auxiliaries")
        if len(arbitrary) != 1:
            raise ValueError(
                f"{name} should contain exactly one ArbBField; found {len(arbitrary)}"
            )
        if mapper_runtime_path is not None and arbitrary[0] != mapper_runtime_path:
            raise ValueError(
                f"{name} ArbBField path {arbitrary[0]!r} != {mapper_runtime_path!r}"
            )
        details[name] = arbitrary[0]

    return {"valid": True, "volumes": details, "count": len(details)}


def ensure_dunendggd(
    destination: os.PathLike | str,
    *,
    install: bool = False,
) -> Path:
    destination = Path(destination).expanduser().resolve()
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", DUNENDGGD_REPOSITORY, str(destination)],
            check=True,
        )

    subprocess.run(["git", "fetch", "--all", "--tags"], cwd=destination, check=True)
    subprocess.run(
        ["git", "checkout", "--detach", DUNENDGGD_COMMIT],
        cwd=destination,
        check=True,
    )

    if install:
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-e", str(destination)],
            check=True,
        )
    return destination


def _target_output_name(
    target: str,
    *,
    tms_shift_mm: float = 10000.0,
    lar_shift_mm: float = 10000.0,
) -> str:
    if target in _TARGET_OUTPUTS:
        return _TARGET_OUTPUTS[target]
    tms = f"{tms_shift_mm:g}"
    lar = f"{lar_shift_mm:g}"
    if target == "prism_nosand":
        return f"nd_hall_with_lar_shift_{lar}_tms_shift_{tms}_nosand.gdml"
    if target == "prism":
        return f"nd_hall_with_lar_shift_{lar}_tms_shift_{tms}_sand_stt1.gdml"
    if target == "prism_drift1":
        return f"nd_hall_with_lar_shift_{lar}_tms_shift_{tms}_sand_drift1.gdml"
    raise ValueError(f"Unsupported dunendggd target {target!r}")


def build_upstream_gdml(
    dunendggd_dir: os.PathLike | str,
    target: str,
    output_path: os.PathLike | str,
    *,
    tms_shift_mm: float = 10000.0,
    lar_shift_mm: float = 10000.0,
) -> Path:
    dunendggd_dir = Path(dunendggd_dir).resolve()
    cmd = ["make", target]
    if target.startswith("prism"):
        cmd += [
            f"TMS_SHIFT={tms_shift_mm:g}",
            f"LAr_SHIFT={lar_shift_mm:g}",
        ]
    subprocess.run(cmd, cwd=dunendggd_dir, check=True)

    generated = dunendggd_dir / _target_output_name(
        target,
        tms_shift_mm=tms_shift_mm,
        lar_shift_mm=lar_shift_mm,
    )
    if not generated.exists():
        raise FileNotFoundError(
            f"dunendggd target {target!r} completed but {generated.name} was not found"
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(generated, output_path)
    return output_path
