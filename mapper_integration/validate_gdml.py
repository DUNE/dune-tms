
#!/usr/bin/env python3
from __future__ import annotations

import argparse

from geometry.gdml_tools import (
    compute_field_alignment,
    inspect_tms_geometry,
    validate_injected_gdml,
)


def main() -> int:
    p = argparse.ArgumentParser(description="Inspect TMS GDML and ArbBField integration.")
    p.add_argument("gdml")
    p.add_argument("--mapper", default=None)
    p.add_argument("--field-dir", default="Field_maps")
    p.add_argument("--mapper-runtime-path", default=None)
    p.add_argument("--tms-volume", default="volTMS")
    p.add_argument("--tolerance-mm", type=float, default=5.0)
    p.add_argument("--source-length-unit", choices=("m", "mm"), default="m")
    args = p.parse_args()

    geom = inspect_tms_geometry(args.gdml, tms_volume=args.tms_volume)
    print("TMS GDML geometry")
    print(f"  global origin [mm] : {geom.global_translation_mm}")
    print(f"  axis aligned       : {geom.axis_aligned}")
    for thickness in (15.0, 40.0, 80.0):
        centers = geom.layer_centers_mm.get(thickness, [])
        print(f"  {thickness:g} mm layers      : {len(centers)}")

    if args.mapper:
        _, report = compute_field_alignment(
            args.gdml,
            args.mapper,
            args.field_dir,
            tms_volume=args.tms_volume,
            source_length_unit=args.source_length_unit,
            tolerance_mm=args.tolerance_mm,
        )
        print("\nField/GDML compatibility")
        print(f"  compatible         : {report.geometry_compatible}")
        print(f"  local translation  : {report.field_to_tms_local_mm} mm")
        print(f"  global translation : {report.field_to_global_mm} mm")
        print(f"  max Z residual     : {report.max_z_residual_mm:.3f} mm")
        for message in report.messages:
            print(f"  - {message}")

    if args.mapper_runtime_path is not None:
        result = validate_injected_gdml(
            args.gdml,
            mapper_runtime_path=args.mapper_runtime_path,
        )
        print("\nArbBField injection")
        print(f"  valid              : {result['valid']}")
        print(f"  volumes            : {result['count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
