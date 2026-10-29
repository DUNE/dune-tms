"""Utilities for building and inspecting TMS magnetic-field mapper files."""

from .core import (
    TMS_COORDINATE_CONVENTION,
    build_uniform_mapper,
    discover_plate_groups,
    read_fld,
    read_z_coordinates,
    reflect_axial_from_quarter,
    source_to_tms_axes,
)
from .edep import read_mapper_header, validate_mapper_file

__all__ = [
    "TMS_COORDINATE_CONVENTION",
    "build_uniform_mapper",
    "discover_plate_groups",
    "read_fld",
    "read_z_coordinates",
    "reflect_axial_from_quarter",
    "source_to_tms_axes",
    "read_mapper_header",
    "validate_mapper_file",
]
