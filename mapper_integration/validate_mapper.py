#!/usr/bin/env python3
from __future__ import annotations

import argparse

from tms_mapper.edep import validate_mapper_file


def main() -> int:
    p = argparse.ArgumentParser(
        description="Validate an ArbBField Mapper.txt for ClarkMcGrew/edep-sim."
    )
    p.add_argument("mapper")
    args = p.parse_args()

    result = validate_mapper_file(args.mapper)
    print("VALID edep-sim ArbBField map")
    print(f"  rows          : {result['rows']:,}")
    print(f"  shape         : {result['shape']}")
    print(f"  offset [mm]   : {result['offset_mm']}")
    print(f"  spacing [mm]  : {result['spacing_mm']}")
    print(f"  max |Bmag-error| [T]: {result['max_bmag_error_t']:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
